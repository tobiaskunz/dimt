/*
 * Copyright (c) 2014, Georgia Tech Research Corporation
 * All rights reserved.
 *
 * Author: Tobias Kunz <tobias@gatech.edu>
 *
 * This file is provided under the following "BSD-style" License:
 *   Redistribution and use in source and binary forms, with or
 *   without modification, are permitted provided that the following
 *   conditions are met:
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 *   CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 *   INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 *   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 *   USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *   AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *   LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *   POSSIBILITY OF SUCH DAMAGE.
 *
 * This file implements a steering method for the double-integrator minimum
 * time problem as described in section IV of [1].
 *
 * [1] Tobias Kunz, Mike Stilman. Probabilistically Complete Kinodynamic
 *     Planning for Robot Manipulators with Acceleration Limits. IEEE/RSJ
 *     International Conference on Intelligent Robots and Systems (IROS), 2014.
 */

#pragma once

#include <iostream>
#include <utility>
#include <array>
#include <cmath>
#include <Eigen/Core>

template <int dof>
class DoubleIntegrator {
public:
	typedef Eigen::Matrix<double, 2 * dof, 1> StateVector;
	typedef Eigen::Matrix<double, dof, 1> Vector;

	class Trajectory1D {
		// 3 segments of constant acceleration: position/velocity[i] := position/velocity at time[i]
		double times[4];
		double positions[4];
		double velocities[4];

	public:

		int getSegment(double t) const {
			assert(0.0 <= t);

			int i = 0;
			while(i < 3 && t >= times[i+1])
				i++;
			return i;
		}

		void getState(double t, double& position, double& velocity) const {
			int i = getSegment(t);
			if(i == 3) {
				position = positions[3];
				velocity = velocities[3];
			}
			else {
				const double s = (t - times[i]) / (times[i+1] - times[i]);
				velocity = (1.0 - s) * velocities[i] + s * velocities[i+1];
				position = (1.0 - s*s) * positions[i] + s*s * positions[i+1] + (t - times[i]) * (1.0 - s) * velocities[i];
			}
			assert(velocity == velocity);
			assert(position == position);
		}
		
		void setStartAndGoal(double startPosition, double startVelocity, double goalPosition, double goalVelocity) {
			positions[0] = startPosition;
			velocities[0] = startVelocity;
			positions[3] = goalPosition;
			velocities[3] = goalVelocity;
		}

		void setPPTrajectory(double duration1, double duration2, double acceleration1) {
			setPLPTrajectory(duration1, 0.0, duration2, velocities[0] + duration1 * acceleration1);
		}

		void setPLPTrajectory(double duration1, double duration2, double duration3, double limitVelocity) {
			// Check for NaNs
			assert(duration1 == duration1);
			assert(duration2 == duration2);
			assert(duration3 == duration3);
			assert(limitVelocity == limitVelocity);

			times[0] = 0.0;
			times[1] = duration1;
			times[2] = duration1 + duration2;
			times[3] = duration1 + duration2 + duration3;

			velocities[1] = velocities[2] = limitVelocity;
			positions[1] = positions[0] + 0.5 * (velocities[0] + velocities[1]) * duration1;
			positions[2] = positions[3] - 0.5 * (velocities[2] + velocities[3]) * duration3;
		}

		double getDuration() const {
			return times[3];
		}

		bool satisfiesVelocityLimits(double minVelocity, double maxVelocity) const {
			return (minVelocity <= velocities[1] && velocities[1] <= maxVelocity);
		}

		bool satisfiesPositionLimits(double minPosition, double maxPosition) const {
			if(velocities[0] * velocities[1] < 0.0) {
				const double startExtremalPosition = positions[0] + 0.5 * velocities[0] * velocities[0] / (velocities[0] - velocities[1]) * times[1];
				if(startExtremalPosition < minPosition || maxPosition < startExtremalPosition)
					return false;
			}

			if(velocities[2] * velocities[3] < 0.0) {
				const double goalExtremalPosition = positions[3] - 0.5 * velocities[3] * velocities[3] / (velocities[3] - velocities[2]) * (times[3] - times[2]);
				if(goalExtremalPosition < minPosition || maxPosition < goalExtremalPosition)
					return false;
			}

			return true;
		}
	};

	struct Trajectory : public std::array<Trajectory1D, dof> {
		double getDuration() const {
			return (*this)[0].getDuration();
		}

		StateVector getState(double time) const {
			StateVector state;
			for(int i = 0; i < dof; i++) {
				(*this)[i].getState(time, state[i], state[dof+i]);
			}
			return state;
		}
	};

private:

	static inline double sign(double d) {
		return d >= 0.0 ? 1.0 : -1.0;
	}

	static inline double squared(double d) {
		return d * d;
	}

	static Trajectory1D getMinAcceleration(double startPosition, double startVelocity, double goalPosition, double goalVelocity, double time, double maxVelocity)
	{
		Trajectory1D trajectory;
		trajectory.setStartAndGoal(startPosition, startVelocity, goalPosition, goalVelocity);
		
		const double a = time * time;
		const double b = 2.0 * (startVelocity + goalVelocity) * time - 4.0 * (goalPosition - startPosition);
		const double c = -squared(goalVelocity - startVelocity);

		if(b == 0.0 && c == 0.0) {
			trajectory.setPLPTrajectory(0.0, time, 0.0, startVelocity);
			return trajectory;
		}

		// numerically stable solution to the quadratic equation
		const double acceleration = (-b - sign(b) * sqrt(b*b - 4.0*a*c)) / (2.0 * a);

		const double t1 = 0.5 * ((goalVelocity - startVelocity) / acceleration + time);
		assert(t1 >= -0.000001);
		assert(t1 <= time + 0.000001);

		trajectory.setPPTrajectory(t1, time - t1, acceleration);

		if(!trajectory.satisfiesVelocityLimits(-maxVelocity, maxVelocity)) {
			// Calculate PLP solution
			const double velocityLimit = sign(acceleration) * maxVelocity;
			if(startVelocity == velocityLimit && goalVelocity == velocityLimit) {
				trajectory.setPLPTrajectory(0.0, time, 0.0, velocityLimit);
				return trajectory;
			}
			const double acceleration = 0.5 * (squared(velocityLimit - startVelocity) + squared(goalVelocity - velocityLimit))
			                            / (velocityLimit * time - (goalPosition - startPosition));
			const double t1 = (velocityLimit - startVelocity) / acceleration;
			const double t2 = (goalVelocity - velocityLimit) / -acceleration;
			assert(t1 == t1);
			trajectory.setPLPTrajectory(t1, time - t1 - t2, t2, velocityLimit);
		}
		return trajectory;
	}

public:
	static double getMinTime1D(double velocity1, double velocity2, double distance,
	                           double maxAcceleration, double maxVelocity, double maxTime,
	                           std::pair<double, double>& infeasibleInterval)
	{
		if(velocity1 > velocity2)
			std::swap(velocity1, velocity2);
		double velocity1Squared = velocity1 * velocity1;
		double velocity2Squared = velocity2 * velocity2;

		const double accelerationDistance = 0.5 * (velocity2Squared - velocity1Squared) / maxAcceleration;
		const double additionalDistance = distance - accelerationDistance;
		
		if(additionalDistance < 0.0)
			std::swap(velocity1Squared, velocity2Squared);

		const double absVelocity = sqrt(maxAcceleration * std::abs(additionalDistance) + velocity2Squared);
		const double velocity = sign(additionalDistance) * absVelocity;
		
		double time = std::abs(velocity1 + velocity2 - 2.0 * velocity) / maxAcceleration;

		if(absVelocity > maxVelocity)
			time += squared(absVelocity - maxVelocity) / (maxVelocity * maxAcceleration);

		if(time >= maxTime)
			return time;

		const double zeroDistance = 0.5 * (velocity1Squared + velocity2Squared) / maxAcceleration;
		if(velocity1 * velocity2 <= 0.0 || additionalDistance * velocity1 < 0.0 || zeroDistance < std::abs(distance)) {
			infeasibleInterval.first = std::numeric_limits<double>::infinity();
			infeasibleInterval.second = std::numeric_limits<double>::infinity();
		}
		else {
			const double absVelocity = sqrt(maxAcceleration * abs(additionalDistance) + velocity1Squared);
			const double velocity = sign(additionalDistance) * absVelocity;
			infeasibleInterval.first = std::abs(velocity1 + velocity2 - 2.0 * velocity) / maxAcceleration;
			infeasibleInterval.second = std::abs(velocity1 + velocity2 + 2.0 * velocity) / maxAcceleration;
		}

		return time;
	}

	static double getMinTime(const StateVector &state1, const StateVector &state2,
	                         const Vector &maxAccelerations,
	                         const Vector &maxVelocities =  std::numeric_limits<double>::infinity() * Vector::Ones(),
	                         double maxTime = std::numeric_limits<double>::infinity())
	{
		const Vector distances = state2.template head<dof>() - state1.template head<dof>();
		double minTime = 0.0;
		std::pair<double, double> infeasibleIntervals[dof];
		
		for(unsigned int i = 0; i < dof && minTime < maxTime; ++i) {
			minTime = std::max(minTime, getMinTime1D(state1[dof+i], state2[dof+i], distances[i], maxAccelerations[i], maxVelocities[i], maxTime, infeasibleIntervals[i]));
			for(unsigned int j = 0; j <= i && minTime < maxTime; ++j) {
				if(infeasibleIntervals[j].first < minTime && minTime < infeasibleIntervals[j].second) {
					minTime = infeasibleIntervals[j].second;
					j = 0;
				}
			}
		}

		assert(minTime < std::numeric_limits<double>::infinity());
		return minTime;
	}

	static Trajectory getTrajectory(const StateVector &state1, const StateVector &state2, double time, const Vector &maxVelocities)
	{
		Trajectory trajectory;
		for(unsigned int i = 0; i < dof; i++)
			trajectory[i] = getMinAcceleration(state1[i], state1[dof + i], state2[i], state2[dof + i], time, maxVelocities[i]);
		return trajectory;
	}

	static Trajectory getTrajectory(const StateVector &state1, const StateVector &state2, const Vector &maxAccelerations, const Vector &maxVelocities)
	{
		const double time = getMinTime(state1, state2, maxAccelerations, maxVelocities);
		return getTrajectory(state1, state2, time, maxVelocities);
	}
};
