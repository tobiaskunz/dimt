#pragma once

#include <Eigen/Core>
#include <iostream>
#include <vector>
#include <fstream>
#include <utility>
#include <Eigen/StdVector>
#include <array>

template <int dof>
class DoubleIntegrator {
public:
	typedef Eigen::Matrix<double, 2 * dof, 1> StateVector;
	typedef Eigen::Matrix<double, dof, 1> Vector;
	typedef Eigen::Matrix<double, dof, dof> Matrix;

	struct Trajectory1D {
		double time[5];
		double limitVelocity;
		double startExtremalPosition;
		double goalExtremalPosition;
		double startPosition;
		double startVelocity;
		double goalPosition;
		double goalVelocity;
		
		int getSegment(double& t) const {
			int i = 0;
			while(t > time[i] && i < 5) {
				t -= time[i];
				i++;
			}
			
			if(i == 5) {
				do {
					i--;
				} while(time[i] == 0.0);
				t = time[i];
			}
			
			assert(t <= time[i]);
			return i;
		}

		void getState(double t, double& position, double& velocity) const {
			int i = getSegment(t);
			assert(t >= 0.0);
			assert(t <= time[i]);
			switch(i) {
			case 0: {
				const double s = t == 0.0 ? 1.0 : (time[0] - t) / time[0];
				velocity = s * startVelocity;
				position = s*s * startPosition + (1 - s*s) * startExtremalPosition;
				break; }
			case 1: {
				const double extremalVelocity = time[0] == 0.0 ? startVelocity : 0.0;
				const double acceleration = (limitVelocity - extremalVelocity) / time[1];
				velocity = extremalVelocity + t * acceleration;
				position = startExtremalPosition + t * extremalVelocity + 0.5 * squared(t) * acceleration;
				break; }
			case 2: {
				const double extremalVelocity = time[0] == 0.0 ? startVelocity : 0.0;
				velocity = limitVelocity;
				position = startExtremalPosition + 0.5 * (limitVelocity + extremalVelocity) * time[1] + limitVelocity * t;
				break; }
			case 3: {
				const double extremalVelocity = time[4] == 0.0 ? goalVelocity : 0.0;
				const double acceleration = (extremalVelocity - limitVelocity) / time[3];
				velocity = extremalVelocity + (t - time[3]) * acceleration;
				position = goalExtremalPosition + (t - time[3]) * extremalVelocity + 0.5 * squared(t - time[3]) * acceleration;
				break; }
			case 4: {
				const double s = t / time[4];
				velocity = s * goalVelocity;
				position = s*s * goalPosition + (1.0 - s*s) * goalExtremalPosition;
				break; }
			}
			assert(position == position);
			assert(velocity == velocity);
		}

		void setStartAndGoal(double startPosition, double startVelocity, double goalPosition, double goalVelocity) {
			this->startPosition = startPosition;
			this->startVelocity = startVelocity;
			this->goalPosition = goalPosition;
			this->goalVelocity = goalVelocity;
		}

		void setPPTrajectory(double time1, double time2, double acceleration1) {
			setPLPTrajectory(time1, 0.0, time2, startVelocity + time1 * acceleration1);
		}

		void setPLPTrajectory(double time1, double time2, double time3, double limitVelocity) {
			// Check for NaNs
			assert(time1 == time1);
			assert(time2 == time2);
			assert(time3 == time3);
			assert(limitVelocity == limitVelocity);

			this->time[2] = time2;
			this->limitVelocity = limitVelocity;

			double distance = goalPosition - startPosition;
			if(startVelocity * limitVelocity >= 0.0) {
				this->startExtremalPosition = startPosition;
				this->time[0] = 0.0;
				this->time[1] = time1;
			}
			else {
				const double stoppingTime = -time1 / (this->limitVelocity - startVelocity) * startVelocity;
				assert(stoppingTime >= 0.0);
				this->startExtremalPosition = startPosition + 0.5 * startVelocity * stoppingTime;
				this->time[0] = stoppingTime;
				this->time[1] = time1 - stoppingTime;
			}

			if(goalVelocity * limitVelocity >= 0.0) {
				this->goalExtremalPosition = goalPosition;
				this->time[3] = time3;
				this->time[4] = 0.0;
			}
			else {
				const double stoppingTime = -time3 / (this->limitVelocity - goalVelocity) * goalVelocity;
				assert(stoppingTime >= 0.0);
				this->goalExtremalPosition = goalPosition - 0.5 * goalVelocity * stoppingTime;
				this->time[3] = time3 - stoppingTime;
				this->time[4] = stoppingTime;
			}
		}

		double getDuration() const {
			return time[0] + time[1] + time[2] + time[3] + time[4];
		}

		bool satisfiesVelocityLimits(double minVelocity, double maxVelocity) const {
			return (minVelocity <= limitVelocity && limitVelocity <= maxVelocity);
		}

		bool satisfiesPositionLimits(double minPosition, double maxPosition) const {
			return (minPosition <= startExtremalPosition && startExtremalPosition <= maxPosition
			     && minPosition <=  goalExtremalPosition &&  goalExtremalPosition <= maxPosition);
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

	static double getPLPTime(double startVelocity, double goalVelocity, double distance, double maxVelocity, double acceleration1, Trajectory1D* trajectory)
	{
		const double boundaryVelocity = sign(acceleration1) * maxVelocity;
		const double timeToBoundary1 = (boundaryVelocity - startVelocity) / acceleration1;
		const double timeToBoundary2 = (goalVelocity - boundaryVelocity) / -acceleration1;
		const double distanceToBoundary = (2.0 * squared(boundaryVelocity) - squared(startVelocity) - squared(goalVelocity)) / 2.0 / acceleration1;
		const double boundaryTime = (distance - distanceToBoundary) / boundaryVelocity;

		if(boundaryTime < 0.0)
			return std::numeric_limits<double>::infinity();
		else {
			if(trajectory)
				trajectory->setPLPTrajectory(timeToBoundary1, boundaryTime, timeToBoundary2, boundaryVelocity);
			return timeToBoundary1 + boundaryTime + timeToBoundary2;
		}
	}

	static void outputTrajectory(std::ostream* output, Trajectory1D* trajectory) {
		//if(trajectory) {
		//	double switch0 = trajectory->time[0];
		//	double switch1 = switch0 + trajectory->time[1];
		//	double switch2 = switch1 + trajectory->time[2];
		//	*output << 0.0 << "  " << trajectory->startState[0] << "  " << trajectory->startState[1] << std::endl
		//		<< switch0 << "  " << trajectory->getState(switch0)[0] << "  " << trajectory->getState(switch0)[1] << std::endl
		//		<< switch1 << "  " << trajectory->getState(switch1)[0] << "  " << trajectory->getState(switch1)[1] << std::endl
		//		<< switch2 << "  " << trajectory->getState(switch2)[0] << "  " << trajectory->getState(switch2)[1] << std::endl
		//		<< switch2 << "  " << trajectory->goalState[0] << "  " << trajectory->goalState[1] << std::endl << std::endl << std::endl;
		//}
		//else {
		//	*output << "0 0 0" << std::endl << std::endl << std::endl;
		//}
	}


	static bool getMinAcceleration(double startPosition, double startVelocity, double goalPosition, double goalVelocity, double time,
	                               double maxAcceleration, double maxVelocity,
	                               Trajectory1D &trajectory)
	{
		trajectory.setStartAndGoal(startPosition, startVelocity, goalPosition, goalVelocity);
		
		const double a = time * time;
		const double b = 2.0 * (startVelocity + goalVelocity) * time - 4.0 * (goalPosition - startPosition);
		const double c = -squared(goalVelocity - startVelocity);

		if(b == 0.0 && c == 0.0) {
			trajectory.setPLPTrajectory(0.0, time, 0.0, startVelocity);
			return true;
		}

		const double radicand = b*b-4.0*a*c;
		assert(radicand >= 0.0);

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
				return true;
			}
			const double acceleration = 0.5 * (squared(velocityLimit - startVelocity) + squared(goalVelocity - velocityLimit))
			                            / (velocityLimit * time - (goalPosition - startPosition));
			const double t1 = (velocityLimit - startVelocity) / acceleration;
			const double t2 = (goalVelocity - velocityLimit) / -acceleration;
			assert(t1 == t1);
			trajectory.setPLPTrajectory(t1, time - t1 - t2, t2, velocityLimit);
		}
		return true;
	}

public:
	
	static std::pair<double, double> getInfeasibleInterval(double startVelocity, double goalVelocity, double distance,
		double acceleration1, double maxAcceleration, double maxVelocity)
	{
		std::pair<double, double> impossibleInterval;

		Trajectory1D boundaryTrajectoryLow;
		boundaryTrajectoryLow.setStartAndGoal(0.0, startVelocity, distance, goalVelocity);
		Trajectory1D boundaryTrajectoryHigh = boundaryTrajectoryLow;
		
		// for all cases below, zero is a feasible acceleration
		// for all cases below: sign(acceleration1) == -sign(acceleration2)
		if(startVelocity * goalVelocity <= 0.0 || acceleration1 * startVelocity < 0.0) {
			// If minimum-time solution goes through zero-velocity, there is no impossible time interval, because we can stop and wait at zero-velocity
			impossibleInterval.first = std::numeric_limits<double>::infinity();
			impossibleInterval.second = std::numeric_limits<double>::infinity();
		}
		// for all cases below: sign(startVelocity) == sign(goalVelocity) == sign(distance) == sign(additionalDistance) == sign(acceleration1) == -sign(acceleration2)
		else {
			double zeroTime1 = abs(startVelocity) / maxAcceleration;
			double zeroTime2 = abs(goalVelocity) / maxAcceleration;
			double zeroDistance = zeroTime1 * startVelocity / 2.0 + zeroTime2 * goalVelocity / 2.0;
			if(abs(zeroDistance) < abs(distance)) {
				impossibleInterval.first = std::numeric_limits<double>::infinity();
				impossibleInterval.second = std::numeric_limits<double>::infinity();
			}
			else {
				double timeLow1, timeLow2, timeHigh1, timeHigh2;
				getPPTime(startVelocity, goalVelocity, distance, -acceleration1, &timeLow1, &timeLow2, &timeHigh1, &timeHigh2);
				impossibleInterval.first = timeLow1 + timeLow2;
				impossibleInterval.second = timeHigh1 + timeHigh2;
				if(impossibleInterval.second == std::numeric_limits<double>::infinity()) {
					std::cout << "infinity 1" << std::endl;
					impossibleInterval.first = 0.0;
					impossibleInterval.second = numeric_limits<double>::infinity();
					return impossibleInterval;
					assert(false);
				}

				if(abs(startVelocity + timeHigh1 * -acceleration1) >= maxVelocity) {
					impossibleInterval.second = getPLPTime(startVelocity, goalVelocity, distance, maxVelocity, -acceleration1, &boundaryTrajectoryHigh);
					if(impossibleInterval.second == std::numeric_limits<double>::infinity()) {
						std::cout << "infinity 2" << std::endl;
						//cout << startVelocity << " " << goalVelocity << " " << distance << " " << minVelocity << " " << maxVelocity << " " << acceleration2 << " " << acceleration1 << endl;
					}
				}
			}
		}

		//if(output) {
		//	if(impossibleInterval.first == std::numeric_limits<double>::infinity())
		//		outputTrajectory(output, NULL);
		//	else
		//		outputTrajectory(output, &boundaryTrajectoryLow);
		//	if(impossibleInterval.second == std::numeric_limits<double>::infinity())
		//		outputTrajectory(output, NULL);
		//	else
		//		outputTrajectory(output, &boundaryTrajectoryHigh);
		//}

		return impossibleInterval;
	}


	static void getPPTime(double startVelocity, double goalVelocity, double distance, double acceleration1,
	                      double* tMin1, double* tMin2, double* tMax1, double* tMax2)
	{
		if(startVelocity == goalVelocity && distance == 0.0) {
			if(tMin1) {
				*tMin1 = 0.0;
				*tMin2 = 0.0;
			}
			if(tMax1) {
				*tMax1 = 0.0;
				*tMax2 = 0.0;
			}
			return;
		}

		const double a = acceleration1;
		const double b = 2.0 * startVelocity;
		//const double c = (squared(goalVelocity) - squared(startVelocity)) / 2.0 / acceleration2 - distance;
		const double c = 0.5 * (startVelocity + goalVelocity) * ((goalVelocity - startVelocity) / -acceleration1) - distance;

		const double radicand = b*b-4.0*a*c;
		assert(radicand >= 0.0);
		
		// numerically stable solution to the quadratic equation
		const double q = -0.5 * (b + sign(b) * sqrt(radicand));

		if(a*sign(b) >= 0.0) {
			if(tMin1) *tMin1 = q / a;
			if(tMax1) *tMax1 = c / q;
		}
		else {
			if(tMin1) *tMin1 = c / q;
			if(tMax1) *tMax1 = q / a;
		}

		if(tMin1) {
			*tMin2 = *tMin1 - (goalVelocity - startVelocity) / acceleration1;
			if(*tMin1 < 0.0 || *tMin2 < 0.0) {
				*tMin1 = numeric_limits<double>::infinity();
				*tMin2 = numeric_limits<double>::infinity();
			}
		}

		if(tMax1) {
			*tMax2 = *tMax1 - (goalVelocity - startVelocity) / acceleration1;
			if(*tMax1 < 0.0 || *tMax2 < 0.0) {
				*tMax1 = numeric_limits<double>::infinity();
				*tMax2 = numeric_limits<double>::infinity();
			}
		}
	}

	static double getMinTime1D(double startVelocity, double goalVelocity, double distance,
	                           double maxAcceleration, double maxVelocity, double maxTime,
	                           double& acceleration1)
	{
		// determine distance travelled while accelerating from start velocity to goal velocity
		const double accelerationTime = abs(goalVelocity - startVelocity) / maxAcceleration;
		if(accelerationTime >= maxTime)
			return accelerationTime;
		const double accelerationDistance = 0.5 * (startVelocity + goalVelocity) * accelerationTime;
		
		// determine whether to accelerate at the upper or lower limit first
		const double additionalDistance = distance - accelerationDistance;

		acceleration1 = sign(additionalDistance) * maxAcceleration;

		double time1, time2;
		getPPTime(startVelocity, goalVelocity, distance, acceleration1, NULL, NULL, &time1, &time2);
		double time = time1 + time2;

		assert(time >= 0.0);
		//assert(time != numeric_limits<double>::infinity());

		if(abs(startVelocity + acceleration1 * time1) >= maxVelocity) {
			//cout << "Hit velocity limit" << endl;
			time = getPLPTime(startVelocity, goalVelocity, distance, maxVelocity, acceleration1, NULL);
		}
		//outputTrajectory(&file, &trajectory);

		return time;
	}


	static void adjustForInfeasibleIntervals(int i, double time, double maxTime,
	                                         const Eigen::Ref<const Vector>& startVelocities, const Eigen::Ref<const Vector>& goalVelocities, const Eigen::Ref<const Vector>& distances,
	                                         const Eigen::Ref<const Vector>& firstAccelerations, const Eigen::Ref<const Vector>& maxAccelerations, const Eigen::Ref<const Vector>& maxVelocities,
	                                         //const Vector& startVelocities, const Vector& goalVelocities, const Vector& distances,
	                                         //const Vector& firstAccelerations, const Vector& maxAccelerations, const Vector& maxVelocities,
	                                         double& minTime, pair<double, double>* infeasibleIntervals, int& limitDof)
	{
		if(time >= maxTime) {
			minTime = time;
			return;
		}

		bool minTimeIncreased = false;

		if(time < minTime) {
			infeasibleIntervals[i] = getInfeasibleInterval(startVelocities[i], goalVelocities[i], distances[i],
			                                               firstAccelerations[i], maxAccelerations[i], maxVelocities[i]);
			if(infeasibleIntervals[i].first < minTime && minTime < infeasibleIntervals[i].second) {
				minTime = infeasibleIntervals[i].second;
				minTimeIncreased = true;
				if(limitDof != -1) {
					infeasibleIntervals[limitDof] = getInfeasibleInterval(startVelocities[limitDof], goalVelocities[limitDof], distances[limitDof],
					                                                      firstAccelerations[limitDof], maxAccelerations[limitDof], maxVelocities[limitDof]);
					limitDof = -1;
				}
			}
		}
		else {
			minTime = time;
			minTimeIncreased = true;
			if(limitDof != -1) {
				infeasibleIntervals[limitDof] = getInfeasibleInterval(startVelocities[limitDof], goalVelocities[limitDof], distances[limitDof],
				                                                      firstAccelerations[limitDof], maxAccelerations[limitDof], maxVelocities[limitDof]);
			}
			limitDof = i;
		}

		while(minTimeIncreased && minTime < maxTime) {
			assert(limitDof == -1 || limitDof == i);
			minTimeIncreased = false;
			for(unsigned int j = 0; j < (limitDof == -1 ? i + 1 : i); j++) {
				if(minTime > infeasibleIntervals[j].first && minTime < infeasibleIntervals[j].second) {
					minTime = infeasibleIntervals[j].second;
					minTimeIncreased = true;
				}
			}
			
			if(limitDof != -1 && minTimeIncreased) {
				infeasibleIntervals[limitDof] = getInfeasibleInterval(startVelocities[limitDof], goalVelocities[limitDof], distances[limitDof],
				                                                      firstAccelerations[limitDof], maxAccelerations[limitDof], maxVelocities[limitDof]);
				limitDof = -1;
			}
		}

		//cout << i << endl;
		//for(int j = 0; j < dof; ++j)
		//	cout << infeasibleIntervals[j].first << "  " << infeasibleIntervals[j].second << endl;
		//cout << endl;
	}

	static double getMinTime(const StateVector &state1, const StateVector &state2,
			const Vector &maxAccelerations,
			const Vector &maxVelocities =  std::numeric_limits<double>::infinity() * Vector::Ones(),
			double maxTime = std::numeric_limits<double>::infinity())
	{
		const VectorBlock<const StateVector, dof> startVelocities = state1.tail<dof>();
		const VectorBlock<const StateVector, dof> goalVelocities = state2.tail<dof>();
		const Vector distances = state2.head<dof>() - state1.head<dof>();

		double minTime = 0.0;
		int limitDof = -1; // DOF for which the min time but not the impossible interval has been calculated yet
		std::pair<double, double> infeasibleIntervals[dof];
		Vector firstAccelerations;
		
		for(unsigned int i = 0; i < dof && minTime < maxTime; ++i) {
			const double time = getMinTime1D(startVelocities[i], goalVelocities[i], distances[i], maxAccelerations[i], maxVelocities[i], maxTime, firstAccelerations[i]);
			adjustForInfeasibleIntervals(i, time, maxTime, startVelocities, goalVelocities, distances, firstAccelerations, maxAccelerations, maxVelocities,
			                             minTime, infeasibleIntervals, limitDof);
		}

		return minTime;
	}

	static bool getTrajectory(const StateVector &state1, const StateVector &state2, double time,
	                          const Vector &maxAccelerations, const Vector &maxVelocities,
	                          Trajectory &trajectory)
	{
		for(unsigned int i = 0; i < dof; i++) {
			if(!getMinAcceleration(state1[i], state1[dof + i], state2[i], state2[dof + i], time, maxAccelerations[i], maxVelocities[i], trajectory[i]))
				return false;
		}
		return true;
	}

	static bool getTrajectory(const StateVector &state1, const StateVector &state2,
	                          const Vector &maxAccelerations, const Vector &maxVelocities,
	                          Trajectory &trajectory)
	{
		const double time = getMinTime(state1, state2, maxAccelerations, maxVelocities);
		if(time == numeric_limits<double>::infinity())
			return false;
		return getTrajectory(state1, state2, time, maxAccelerations, maxVelocities, trajectory);
	}
};