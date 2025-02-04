#!/usr/bin/env python3
import math
import numpy as np
from common.numpy_fast import interp

import cereal.messaging as messaging
from common.filter_simple import FirstOrderFilter
from common.realtime import DT_MDL
from selfdrive.modeld.constants import T_IDXS
from selfdrive.config import Conversions as CV
from selfdrive.controls.lib.longcontrol import LongCtrlState
from selfdrive.controls.lib.longitudinal_mpc_lib.long_mpc import LongitudinalMpc
from selfdrive.controls.lib.longitudinal_mpc_lib.long_mpc import T_IDXS as T_IDXS_MPC
from selfdrive.controls.lib.drive_helpers import V_CRUISE_MAX, CONTROL_N
from selfdrive.swaglog import cloudlog
from selfdrive.car.tesla.interface import get_tesla_accel_limits
from selfdrive.car.modules.CFG_module import load_bool_param,load_float_param


LON_MPC_STEP = 0.2  # first step is 0.2s
AWARENESS_DECEL = -0.2  # car smoothly decel at .2m/s^2 when user is distracted
#THESE ARE NOT USED HERE, NEW FUNCTION IN TESLA INTERFACE
#A_CRUISE_MIN = -1.2 
#A_CRUISE_MAX_VALS = [1.2, 1.2, 0.8, 0.6]
#A_CRUISE_MAX_BP = [0., 7.5, 15., 25., 40.]

# Lookup table for turns
_A_TOTAL_MAX_V = [2.2, 4.15]
_A_TOTAL_MAX_BP = [20., 40.]

ACCEL_MIN_TURN_SLOWDOWN = - 1.0 # m/s^2


def get_max_accel(CP,v_ego):
  return get_tesla_accel_limits(CP,v_ego)  

def limit_accel_in_turns(v_ego, angle_steers, a_target, CP):
  """
  This function returns a limited long acceleration allowed, depending on the existing lateral acceleration
  this should avoid accelerating when losing the target in turns
  """

  a_total_max = interp(v_ego, _A_TOTAL_MAX_BP, _A_TOTAL_MAX_V)
  a_y = v_ego ** 2 * angle_steers * CV.DEG_TO_RAD / (CP.steerRatio * CP.wheelbase)
  a_x_allowed = math.sqrt(max(a_total_max ** 2 - a_y ** 2, 0.))

  return [a_target[0], min(a_target[1], a_x_allowed)]


class Planner:
  def __init__(self, CP, init_v=0.0, init_a=0.0):
    self.CP = CP
    self.mpc = LongitudinalMpc()

    self.fcw = False

    self.a_desired = init_a
    self.v_desired_filter = FirstOrderFilter(init_v, 2.0, DT_MDL)

    self.v_desired_trajectory = np.zeros(CONTROL_N)
    self.a_desired_trajectory = np.zeros(CONTROL_N)
    self.j_desired_trajectory = np.zeros(CONTROL_N)

    #used for slow down in turns
    self.leadsData = None
    self.path_x = np.arange(192)
    self.enable_turn_slowdown = load_bool_param("TinklaTurnSlowdown", True)
    self.turn_slowdown_factor = load_float_param("TinklaTurnSlowdownFactor",0.95)

  def get_path_length_idx(self, y, distance):
    i = 0
    for val in y:
        if val < distance:
            i = i + 1
    return i

  def update(self, sm):
    v_ego = sm['carState'].vEgo
    a_ego = sm['carState'].aEgo

    v_cruise_kph = sm['controlsState'].vCruise
    v_cruise_kph = min(v_cruise_kph, V_CRUISE_MAX)
    v_cruise = v_cruise_kph * CV.KPH_TO_MS

    long_control_state = sm['controlsState'].longControlState
    force_slow_decel = sm['controlsState'].forceDecel

    if sm['radarState'] is not None:
      self.leadsData = sm['radarState']

    if (sm['modelV2'] is not None) and self.enable_turn_slowdown:
      #TODO: Use probability to decide if the speed limit should be valid
      #leftLaneQuality = 1 if sm['modelV2'].laneLineProbs[0] > 0.25 else 0
      #rightLaneQuality = 1 if sm['modelV2'].laneLineProbs[3] > 0.25 else 0
      #let's get the position points and compute poly coef
      y = np.array(sm['modelV2'].position.y)
      x = np.array(sm['modelV2'].position.x)
      max_distance = 100.0
      if self.leadsData is not None:
          if self.leadsData.leadOne.status:
              lead_d = self.leadsData.leadOne.dRel * 2.0
              max_distance = max(0,min(lead_d, max_distance))
      max_idx = self.get_path_length_idx(y, max_distance)
      order = 3
      coefs = np.polyfit(x[:max_idx], y[:max_idx], order)
      # Curvature of polynomial https://en.wikipedia.org/wiki/Curvature#Curvature_of_the_graph_of_a_function
      # y = a x^3 + b x^2 + c x + d, y' = 3 a x^2 + 2 b x + c, y'' = 6 a x + 2 b
      # k = y'' / (1 + y'^2)^1.5
      # TODO: compute max speed without using a list of points and without numpy
      y_p = 3 * coefs[0] * self.path_x ** 2 + 2 * coefs[1] * self.path_x + coefs[2]
      y_pp = 6 * coefs[0] * self.path_x + 2 * coefs[1]
      curv = y_pp / (1. + y_p ** 2) ** 1.5
      a_y_max = 3.1 - v_ego * 0.032
      v_curvature = np.sqrt(a_y_max / np.clip(np.abs(curv), 1e-4, None))
      model_speed = np.min(v_curvature) * self.turn_slowdown_factor
      model_speed = max(20.0 * CV.MPH_TO_MS, model_speed)  # Don't slow down below 20mph
    else:
      model_speed = 255.  # (MAX_SPEED)
    #print("curvature_speed=",model_speed, " cruise_speed=",v_cruise )
    #force the speed to the min between what's set and what we need for curvature
    slowdown_for_turn = False
    if model_speed < v_cruise:
      v_cruise = min(v_cruise,model_speed)
      slowdown_for_turn = True

    prev_accel_constraint = True
    if long_control_state == LongCtrlState.off or sm['carState'].gasPressed:
      self.v_desired_filter.x = v_ego
      self.a_desired = a_ego
      # Smoothly changing between accel trajectory is only relevant when OP is driving
      prev_accel_constraint = False

    # Prevent divergence, smooth in current v_ego
    self.v_desired_filter.x = max(0.0, self.v_desired_filter.update(v_ego))

    accel_limits = get_max_accel(self.CP,v_ego)

    accel_limits_turns = limit_accel_in_turns(v_ego, sm['carState'].steeringAngleDeg, accel_limits, self.CP)
    if force_slow_decel:
      # if required so, force a smooth deceleration
      accel_limits_turns[1] = min(accel_limits_turns[1], AWARENESS_DECEL)
      accel_limits_turns[0] = min(accel_limits_turns[0], accel_limits_turns[1])
    # clip limits, cannot init MPC outside of bounds
    accel_limits_turns[0] = min(accel_limits_turns[0], self.a_desired + 0.05)
    accel_limits_turns[1] = max(accel_limits_turns[1], self.a_desired - 0.05)
    self.mpc.set_accel_limits(accel_limits_turns[0], accel_limits_turns[1])
    self.mpc.set_cur_state(self.v_desired_filter.x, self.a_desired)
    self.mpc.update(sm['carState'], sm['radarState'], v_cruise, prev_accel_constraint=prev_accel_constraint)
    self.v_desired_trajectory = np.interp(T_IDXS[:CONTROL_N], T_IDXS_MPC, self.mpc.v_solution)
    self.a_desired_trajectory = np.interp(T_IDXS[:CONTROL_N], T_IDXS_MPC, self.mpc.a_solution)
    self.j_desired_trajectory = np.interp(T_IDXS[:CONTROL_N], T_IDXS_MPC[:-1], self.mpc.j_solution)

    # TODO counter is only needed because radar is glitchy, remove once radar is gone
    self.fcw = self.mpc.crash_cnt > 5
    if self.fcw:
      cloudlog.info("FCW triggered")

    # Interpolate 0.05 seconds and save as starting point for next iteration
    a_prev = self.a_desired
    self.a_desired = float(interp(DT_MDL, T_IDXS[:CONTROL_N], self.a_desired_trajectory))
    if slowdown_for_turn and (v_ego > model_speed):
      #we need to slow down, but not faster than ACCEL_MIN_TURN_SLOWDOWN
      self.a_desired = min(self.a_desired,max (ACCEL_MIN_TURN_SLOWDOWN,model_speed - v_ego))

    self.v_desired_filter.x = self.v_desired_filter.x + DT_MDL * (self.a_desired + a_prev) / 2.0

  def publish(self, sm, pm):
    plan_send = messaging.new_message('longitudinalPlan')

    plan_send.valid = sm.all_alive_and_valid(service_list=['carState', 'controlsState'])

    longitudinalPlan = plan_send.longitudinalPlan
    longitudinalPlan.modelMonoTime = sm.logMonoTime['modelV2']
    longitudinalPlan.processingDelay = (plan_send.logMonoTime / 1e9) - sm.logMonoTime['modelV2']

    longitudinalPlan.speeds = self.v_desired_trajectory.tolist()
    longitudinalPlan.accels = self.a_desired_trajectory.tolist()
    longitudinalPlan.jerks = self.j_desired_trajectory.tolist()

    longitudinalPlan.hasLead = sm['radarState'].leadOne.status
    longitudinalPlan.longitudinalPlanSource = self.mpc.source
    longitudinalPlan.fcw = self.fcw

    longitudinalPlan.solverExecutionTime = self.mpc.solve_time

    pm.send('longitudinalPlan', plan_send)
