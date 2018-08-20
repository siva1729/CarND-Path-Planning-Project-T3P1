
### Introduction ###
---
In this project goal is to safely navigate around a virtual highway with other traffic that is driving +-10 MPH of the 50 MPH speed limit. You will be provided the car's localization and sensor fusion data, there is also a sparse map list of way-points around the highway. The car should try to go as close as possible to the 50 MPH speed limit, which means passing slower traffic when possible, note that other cars will try to change lanes too. The car should avoid hitting other cars at all cost as well as driving inside of the marked road lanes at all times, unless going from one lane to another. The car should be able to make one complete loop around the 6946m highway. Since the car is trying to go 50 MPH, it should take a little over 5 minutes to complete 1 loop. Also the car should not experience total acceleration over 10 m/s^2 and jerk that is greater than 10 m/s^3.

#### The map of the highway is in data/highway_map.txt
Each way-point in the list contains  [x,y,s,dx,dy] values. x and y are the way-point's map coordinate position, the s value is the distance along the road to get to that way-point in meters, the dx and dy values define the unit normal vector pointing outward of the highway loop.

The highway's way-points loop around so the frenet s value, distance along the road, goes from 0 to 6945.554.
The goal of this Path Planning project is to navigate a car around a simulated highway scenario, including traffic and given way-points, telemetry, and sensor fusion data. The car must not violate a set of motion constraints, namely maximum velocity, maximum acceleration, and maximum jerk, while also avoiding collisions with other vehicles, keeping to within a highway lane (aside from short periods of time while changing lanes), and changing lanes when doing so is necessary to maintain a speed near the posted speed limit.

### Implementation ###
---
The path planning project is implemented using the following steps:

1. Generate Predictions from Sensor Fusion Data 
2. Behavior Planning, Generate trajectories using Predictions
3. Determine Best Trajectory using Cost functions
4. Create way-points for new path

![](./img1.jpg)

#### 1. Generate Predictions from Sensor Fusion Data ####
At first, using the Sensor fusion data the behavior and state of the surrounding vehicles/cars are detected (code is presented in line numbers 309 to 398). The stats are calculated for left/right/current lanes during each time data is presented. When the car is in left most lane, the left valid is marked invalid lane and similarly for right most lane. The sensor fusion data typically had captured data for 12 vehicles for which it provided the vehicles x-position, y-position, velocity  components (vx, vy) and frenet co-ordinates s and d.  Using this data, for each vehicles the car's position (lane number), speed (converted to mph)is calculated. The car's position is calculated using the frenet coordinate the lateral displacement (d) which is converted to lane number using the given lane width (4 m). The longitudinal displacement (s), relative to the ego vehicle provided the buffer distance that the ego vehicle has w.r.t to the closet car. The car's speed is calculated using the velocity complements and converted to mph (miles per hour). Using each vehicle's position (lane position and longitudinal displacement) and speed, the following are calculated:

1. average speed of each lane
2. buffer distance between ego vehicle and the closest vehicles in all 3 lanes
3. speed of closest vehicle

			// Check for cars in left/right lanes
			for (int i=0; i < sensor_fusion.size(); i++)
			{
				double vx = sensor_fusion[i][3];
				double vy = sensor_fusion[i][4];
				double check_speed = sqrt(vx*vx+vy*vy);
				double check_car_s = sensor_fusion[i][5];
				
				// Get Lane Number
				int lane_num = getLaneNum(sensor_fusion[i][6]);
				if (lane_num < 0) continue;
								
				// if using previous points can project s value out
				check_car_s +=((double)prev_size*.02*check_speed);
				
				// Get the lane type in relation to current lane
				int check_lane = getLaneType(current_lane, lane_num); //0: same, -1: left, +1: right
				
				// Ignore if check_car belongs to far-off lane
				if (check_lane == NL) continue; 
				
				// Accumulate lane speeds for avg. speed calc.
                check_speed *= 2.24; // conversion to mph
				avgspeed_in_lanes[check_lane] += check_speed; // avg. speed calc.
				count_in_lanes[check_lane] += 1; // count cars.
				
				// Check ahead(front)only in Given Lane				
				if (check_lane == GL) back_buf_dist = 0;
                else                  back_buf_dist = BACK_BUF_DIST;
				
				if (((car_s - back_buf_dist) < check_car_s) && ((car_s + front_buf_dist) > check_car_s))
					cars_in_lanes[check_lane] = true; // car in left lan
				
				// Debug - record only ahead car for GL
				double buf_dist = check_car_s - car_s;
                
                // record the closest vehicle/car either front/back
				if (fabs(buf_dist) < fabs(buffer_in_lanes[check_lane]) && (check_lane != GL || buf_dist > 0)) 
                {
           			buffer_in_lanes[check_lane] = buf_dist;  // check closest
                }
                // record back vehicles/cars
                if (buf_dist < 0 && fabs(buf_dist) < fabs(backbuf_in_lanes[check_lane])) 
                {
                    backbuf_in_lanes[check_lane] = buf_dist;
                    maxspeed_in_lanes[check_lane] = check_speed;
                }
			} 
			

#### 2. Behavior Planning, Generate trajectories using Predictions ####
These predictions are input to the Behavior planning that uses cost functions with appropriate weights to calculate the next best lane suitable for the ego vehicle to take (code is in line numbers 410 to 447) . The cost functions are provided for the ego vehicle to
- maintain the target speed of 49.5 mph (max speed is 50 mph)
- change lanes safely  to navigate/maintain target speed
- to maintain comfort, cost for frequently changing the lanes

####3. Determine Best Trajectory using Cost functions (CF): ####
The cost function (CF\_TARGET\_SPEED) for maintaining the target speed  is calculated based on the average speed for each lane and it changes linearly based on the difference between target speed and average speed of the lane. The CF is normalized to be in the range of [0,1). Similarly a constant cost value (1) is added for changing lanes i.e. this is the cost for moving into left/right lanes (CF\_CHANGE\_LANES). And based on the buffer distance if there is a vehicle within certain range here (+/- 30 m), cost (of value 1) is added for each lane (CF\_BUFFER\_DIST). The weight factors for each cost function are chosen in the following order
CF\_BUFFER\_DIST > CF\_TARGET\_SPEED > CF\_CHANGE\_LANES, so that safety is given the highest priority to avoid collision, then keeping up the target speed along with comfort of not changing lanes frequently.
Also the duration for the ego vehicle is calculated/maintained in a given lane after changing lanes. This was done
to make sure the ego vehicle stays in a lane after changing lanes for certain period of time before changing lanes again. This would also prevent frequent lane changes and maintain safety adding to comfort/safety factors.

            const double weight_lane_change = 20;  // penalty for changing lanes
			const double weight_speed = pow(10,2); // penalty for not keeping avg. speed
			const double weight_safety = pow(10,3); // penalty for safety cocern
			
            const int best_lane_count_limit = 20; // min. # to Stay in given lane before changing lanes

			const double max_acc = 0.224;
			const double target_speed = max_vel;
			double best_cost = numeric_limits<double>::max();
			
			int next_best_lane = GL;
			
			
			for (int lane = LL; lane <= RL; lane++)
			{
				if (count_in_lanes[lane] == -1) continue; // Invalid lane
				double cost = 0; // initial cost
				
				// cost for changing lanes
				if (lane != GL) 
					cost += (1.0 * weight_lane_change);
				
				// cost for not keeping the ref_vel
				cost += (speed_cost(target_speed, avgspeed_in_lanes[lane]) * weight_speed);
									
				// cost for Safety
				if (cars_in_lanes[lane])
					cost += (1.0 * weight_safety);
                
				cost_lanes[lane] = cost;
				
				// Chose the next_lane based on cost
				if (cost < best_cost) 
				{
					next_best_lane = lane;
					best_cost = cost;
				}
			}		

Once the best trajectory is determined for the ego vehicle either to stay in same lane or to change lanes (left/right), then decision is made about changing the velocity. If ego is determined to stay in same lane and there is no vehicle ahead then speed is increased to meet target speed otherwise speed is reduced to avoid collision with vehicle ahead. If ego vehicle is determined to change lanes then lanes are changes when it is safe to do so (code is given with line numbers 450 to 468)

    // Stay in the current/given lane for certain time before changing lanes
			++best_lane_count[next_best_lane];
			if (best_lane_count[next_best_lane] > best_lane_count_limit) 
            {
				if (next_best_lane != GL && cars_in_lanes[next_best_lane] == false) 
                    change_lanes = true;
                //reset best_lane_count: best_lane_count = {0,0,0}; 
                std::fill(best_lane_count.begin(), best_lane_count.end(), 0);
            }
			
			// If too_close and staying in same lane reduce speed else increase speed			
			//if (cars_in_lanes[GL] && (next_best_lane == GL))
            if (cars_in_lanes[GL] && !change_lanes)
			{
				ref_vel -= max_acc;
			}
			else if (ref_vel < max_vel)
			{
				ref_vel += max_acc;
			}

####4. Create way-points for new path ####
 To estimate the location of points between the known way-points, we will need to "interpolate" the position of those points. Fitting polynomials to way-points can be done using quintic polynomial fitting but have chosen spline fitting and which guarantees that the generated function passes through every point. It was easy to use spline function using the code for spline function take from following location [http://kluge.in-chemnitz.de/opensource/spline/](http://kluge.in-chemnitz.de/opensource/spline/)



### Conclusion: ###
---
 The code was able to compile without any errors and the ego vehicle was able to complete more than 1 loop (6564 m) without any issues (i.e. not breaking any requirements) around 5 min and 45 sec. The code made sure the target speed of 49.5 mph was attempted taking into account safety and comfort concerns. It was a very challenging project and was fun after some time when initial requirements were met and code was improved to achieve target speed most of the times without any violations (max speed, jerk, max acceleration,..). For real time scenarios I think human behavior of considering lane change, stability could be considered to improve the code performance.

![](./img2.jpg)
### Basic Build Instructions
---

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./path_planning`.

Here is the data provided from the Simulator to the C++ Program

##### Main car's localization Data (No Noise)

["x"] The car's x position in map coordinates

["y"] The car's y position in map coordinates

["s"] The car's s position in frenet coordinates

["d"] The car's d position in frenet coordinates

["yaw"] The car's yaw angle in the map

["speed"] The car's speed in MPH

##### Previous path data given to the Planner

//Note: Return the previous list but with processed points removed, can be a nice tool to show how far along
the path has processed since last time. 

["previous_path_x"] The previous list of x points previously given to the simulator

["previous_path_y"] The previous list of y points previously given to the simulator

##### Previous path's end s and d values 

["end_path_s"] The previous list's last point's frenet s value

["end_path_d"] The previous list's last point's frenet d value

##### Sensor Fusion Data, a list of all other car's attributes on the same side of the road. (No Noise)

["sensor_fusion"] A 2d vector of cars and then that car's [car's unique ID, car's x position in map coordinates, car's y position in map coordinates, car's x velocity in m/s, car's y velocity in m/s, car's s position in frenet coordinates, car's d position in frenet coordinates. 

#### Details
---

1. The car uses a perfect controller and will visit every (x,y) point it recieves in the list every .02 seconds. The units for the (x,y) points are in meters and the spacing of the points determines the speed of the car. The vector going from a point to the next point in the list dictates the angle of the car. Acceleration both in the tangential and normal directions is measured along with the jerk, the rate of change of total Acceleration. The (x,y) point paths that the planner recieves should not have a total acceleration that goes over 10 m/s^2, also the jerk should not go over 50 m/s^3. (NOTE: As this is BETA, these requirements might change. Also currently jerk is over a .02 second interval, it would probably be better to average total acceleration over 1 second and measure jerk from that.

2. There will be some latency between the simulator running and the path planner returning a path, with optimized code usually its not very long maybe just 1-3 time steps. During this delay the simulator will continue using points that it was last given, because of this its a good idea to store the last points you have used so you can have a smooth transition. previous_path_x, and previous_path_y can be helpful for this transition since they show the last points given to the simulator controller with the processed points already removed. You would either return a path that extends this previous path or make sure to create a new path that has a smooth transition with this last path.

#### Dependencies
---

* cmake >= 3.5
  * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools]((https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
* [uWebSockets](https://github.com/uWebSockets/uWebSockets)
  * Run either `install-mac.sh` or `install-ubuntu.sh`.
  * If you install from source, checkout to commit `e94b6e1`, i.e.
    ```
    git clone https://github.com/uWebSockets/uWebSockets 
    cd uWebSockets
    git checkout e94b6e1
    ```

