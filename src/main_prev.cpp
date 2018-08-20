#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"
#include "math.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
    if (closestWaypoint == maps_x.size())
    {
      closestWaypoint = 0;
    }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s, frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

//========== Helper functions =====================
#define DEBUG_ENABLE

//Lane Type: No Lane, Left Lane, Given Lane, Right Lane
enum lane_type:int {NL = -1, LL = 0, GL, RL};

// constants for lane structures
const int RML = 2; // Right Most Lane number
const int LML = 0; // Left Most Lane number
			
// Cal. the lane number based on Frenet co-ordinates
int getLaneNum (double d) 
{ 	
    const double lane_width = 4;
	int lane = -1;
	if (d > 0 && d < lane_width) 
		lane = 0;
	else if (d > 1*lane_width && d < 2*lane_width) 
		lane = 1;
	else if (d > 2*lane_width && d < 3*lane_width) 
		lane = 2;
	
	return lane;
}

int getLaneType (int current_lane, int given_lane) 
{ 	
    int lane = NL;
	int lane_diff = current_lane - given_lane;
	if      (lane_diff ==  0) lane = GL;
	else if (lane_diff == +1) lane = LL;
	else if (lane_diff == -1) lane = RL;
	
	return lane;
}

double speed_cost(double target_speed, double intended_speed) 
{
	double cost = (1.0*target_speed - intended_speed)/target_speed;
	return cost;
}


int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  // start in lane 1
  int current_lane = 1;
  
  // Reference velocity to target
  double ref_vel = 0; // starting velocity
  double max_vel = 49.5;
  
  
  h.onMessage(
      [&current_lane, &ref_vel, &max_vel, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy]
       (uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
          uWS::OpCode opCode) {
        // "42" at the start of the message means there's a websocket message event.
        // The 4 signifies a websocket message
        // The 2 signifies a websocket event
        //auto sdata = string(data).substr(0, length);
        //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry") {
          // j[1] is the data JSON object
          		
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];
			json msgJson;
			
			// Get previous path size
			int prev_size = previous_path_x.size(); 
			
			

			// Logic to avoid collision vehicles ahead of the Cartesian
			if (prev_size > 0)
			{
				car_s = end_path_s;
			}
			
			// Predictions:
			// Calculate/Predict the Average Lane velocities based on SensorFusion Data for each lane
			// Initialize the lane stats
			
			/* for (int i=0; i < sensor_fusion.size(); i++)
			{
				float d = sensor_fusion[i][6];
				if (d < (2+4*current_lane+2) && d > (2+4*current_lane-2)) 
				{
					double vx = sensor_fusion[i][3];
					double vy = sensor_fusion[i][4];
					double check_speed = sqrt(vx*vx+vy*vy);
					double check_car_s = sensor_fusion[i][5];
					
					// if using previous points can project s value out
					check_car_s +=((double)prev_size*.02*check_speed);
					// check s values greater than mine and s gap
					if ((check_car_s > car_s) && ((check_car_s - car_s) < 30))
					{
						// Flag too close
						too_close = true;
						speed_ahead = check_speed;
						//if (current_lane > 0) { current_lane = 0; }
					}
				}
			} */
			
			vector<bool> cars_in_lanes = {false, false, false}; // {Left, Ahead, Right}
			vector<double> speed_in_lanes = {0, 0, 0}; // {Left, Ahead, Right}
			vector<int> count_in_lanes = {0, 0, 0}; // {Left, Ahead, Right}
			vector<double> buffer_in_lanes = {9999, 9999, 9999}; // {Left, Ahead, Right}
			double back_buf_dist = 30;
			double front_buf_dist = 30;
			
			// Check for invalid lanes if you are in Left-Most/Right-Most lanes
			if (current_lane == RML) 
			{
				cars_in_lanes[RL] = true; 
				count_in_lanes[RL] = -1;
			}
			if (current_lane == LML) 
			{
				cars_in_lanes[LL] = true; 
				count_in_lanes[LL] = -1;
			}			
			
			bool too_close = false;
			
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
				speed_in_lanes[check_lane] += check_speed; // avg. speed calc.
				count_in_lanes[check_lane] += 1; // count cars.
				
				// Check only ahead in Given Lane				
				if (check_lane == GL) back_buf_dist = 0;
				
				if (((car_s - back_buf_dist) < check_car_s) && ((car_s + front_buf_dist) > check_car_s))
					cars_in_lanes[check_lane] = true; // car in left lan
				
				// Debug - record only ahead car for GL
				double buf_dist = check_car_s - car_s;
				if (fabs(buf_dist) < fabs(buffer_in_lanes[check_lane]) && (check_lane != GL || buf_dist > 0))
					buffer_in_lanes[check_lane] = buf_dist;
			} 
			
			// Cal. average speed of lanes {LL, GL, RL}
			for (int lane = LL; lane <= RL; lane++)
			{
				if (count_in_lanes[lane] == -1)     // Invalid lane
					speed_in_lanes[lane] = -1;     // small number
				else if (count_in_lanes[lane] == 0) // No cars in the lane
					speed_in_lanes[lane] = max_vel;
				else 
				{
					speed_in_lanes[lane] *= 2.24; // conversion to mph
					speed_in_lanes[lane] /= ((double)count_in_lanes[lane]); //average speed
				}
			}
			
					
			// CostFunction: Determine the Next Best Lane
			// Calculate associated cost with each lane based on the predictions made so far
			const double max_cost = numeric_limits<double>::max();
			vector<int> cost_lanes = {-1, -1, -1};
			static vector<int> best_lane_count = {0, 0, 0};
			bool change_lanes = false;
			int best_lane_count_limit = 20;
			
			// Weights for the cost calculation
			const double weight_lane_change = 20;
			const double weight_speed = pow(10,2);
			const double weight_safety = pow(10,3);
			
			const double max_acc = 0.224;
			const double target_speed = max_vel;
			double best_cost = numeric_limits<double>::max();
			
			int next_best_lane = GL;
			
			// Keep track of time before changing lanes to avoid switching between lanes abruptly
			const int keep_lane_limit = 200;
			static int keep_lane_count = 0;
			
			for (int lane = LL; lane <= RL; lane++)
			{
				if (count_in_lanes[lane] == -1) continue; // Invalid lane
				double cost = 0; // initial cost
				
				// cost for changing lanes
				if (lane != GL) 
					cost += (1.0 * weight_lane_change);
				
				// cost for not keeping the ref_vel
				cost += (speed_cost(target_speed, speed_in_lanes[lane]) * weight_speed);
									
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
			
			// Keep track of lane changing
			++best_lane_count[next_best_lane];
			if (next_best_lane == RL) best_lane_count[LL] = 0;
			if (next_best_lane == LL) best_lane_count[RL] = 0;
			
			if (best_lane_count[next_best_lane] > best_lane_count_limit)
				change_lanes = true;
			
			if (next_best_lane == GL) ++keep_lane_count;
			else keep_lane_count = 0;
			
			// If too_close and staying in same lane reduce speed else increase speed			
			if (cars_in_lanes[GL] && (next_best_lane == GL))
			{
				ref_vel -= max_acc;
			}
			else if (ref_vel < max_vel)
			{
				ref_vel += max_acc;
			}
			
#ifdef DEBUG_ENABLE			
			//***** Debug ********************************
			int static lane_change_count = 0;
			bool CLL = cars_in_lanes[0];
			bool CGL = cars_in_lanes[1];
			bool CRL = cars_in_lanes[2];
			
			int SLL = speed_in_lanes[0];
			int SGL = speed_in_lanes[1];
			int SRL = speed_in_lanes[2];
			
			// int NLL = count_in_lanes[0];
			// int NGL = count_in_lanes[1];
			// int NRL = count_in_lanes[2];
			
			int NLL = best_lane_count[0];
			int NGL = best_lane_count[1];
			int NRL = best_lane_count[2];

			double PLL = cost_lanes[0];
			double PGL = cost_lanes[1];
			double PRL = cost_lanes[2];
			
			int BLL = buffer_in_lanes[0];
			int BGL = buffer_in_lanes[1];
			int BRL = buffer_in_lanes[2];
			
			if (next_best_lane != GL) lane_change_count++;
			
			int next_lane = current_lane;
			if (next_best_lane == RL) next_lane = current_lane + 1;
			else if (next_best_lane == LL) next_lane = current_lane - 1;
			
			cout << "lane:" << current_lane << "->" << next_lane \
			<<"(" << lane_change_count << ")" << " V: " << int(ref_vel) \
			<< " KLC: " << keep_lane_count \
			<< " Ahead: " << CLL << "," << CGL << "," << CRL \
			<< " Vel: " << SLL << "," << SGL << "," << SRL \
			<< " Buf: " << BLL << "," << BGL << "," << BRL \
			<< " BLC: " << NLL << "," << NGL << "," << NRL \
			<< " cost: " << PLL << "," << PGL << "," << PRL \
			<< endl;
			//************************************************
#endif
			
			//**********************************
			// Change lanes based on cost value:
			// convert to lane number
			if (change_lanes)
			{
				if (next_best_lane == RL) current_lane = current_lane + 1;
				else if (next_best_lane == LL) current_lane = current_lane - 1;
			}
			
			// Create a list of widely spaced (x,y) waypoints, evenly spaced at 30m
			// Then interpolate these waypoints with a spline and fill it with more points that control speed
          	vector<double> ptsx;
          	vector<double> ptsy;
			
			// reference x,y, yaw states
			// use starting point or previous path end point as reference
			double ref_x = car_x;
			double ref_y = car_y;
			double ref_yaw = deg2rad(car_yaw);
			
			//1. Add 2 waypoints to ptsx,ptsy using the previous path way points
			// Check is previous size is almost empty, then use car current position as starting reference
			if (prev_size < 2) 
			{
				//Use two points that make the path tangent to the car
				// Here assumption the distance between points is "1" unit
				double prev_car_x = car_x - cos(car_yaw);
				double prev_car_y = car_y - sin(car_yaw);
				
				ptsx.push_back(prev_car_x);
				ptsx.push_back(car_x);
				ptsy.push_back(prev_car_y);
				ptsy.push_back(car_y);
			}
			else
			{
				// redefine the reference state as previous path end point
				ref_x = previous_path_x[prev_size-1];
				ref_y = previous_path_y[prev_size-1];
				
				double ref_x_prev = previous_path_x[prev_size-2];
				double ref_y_prev = previous_path_y[prev_size-2];
				// Use two pints that mak the path tangent to prev path's end point
				ref_yaw = atan2(ref_y-ref_y_prev, ref_x-ref_x_prev);
				ptsx.push_back(ref_x_prev);
				ptsx.push_back(ref_x);
				ptsy.push_back(ref_y_prev);
				ptsy.push_back(ref_y);
			}
			
			
			//2. Add 3 points at distane of 30,60,90m ahead of starting reference point
            vector<double> next_wp0 = getXY(car_s+30,(2+4*(current_lane)),map_waypoints_s,map_waypoints_x,map_waypoints_y);
			vector<double> next_wp1 = getXY(car_s+60,(2+4*(current_lane)),map_waypoints_s,map_waypoints_x,map_waypoints_y);
			vector<double> next_wp2 = getXY(car_s+90,(2+4*(current_lane)),map_waypoints_s,map_waypoints_x,map_waypoints_y);

			ptsx.push_back(next_wp0[0]);
			ptsx.push_back(next_wp1[0]);
			ptsx.push_back(next_wp2[0]);
			ptsy.push_back(next_wp0[1]);
			ptsy.push_back(next_wp1[1]);
			ptsy.push_back(next_wp2[1]);
			
			// Shift the car reference angle to 0 degrees
			for (int i=0; i < ptsx.size(); i++) 
			{
				double shift_x = ptsx[i]-ref_x;
				double shift_y = ptsy[i]-ref_y;
				
				ptsx[i] = (shift_x*cos(0-ref_yaw)-shift_y*sin(0-ref_yaw));
				ptsy[i] = (shift_x*sin(0-ref_yaw)+shift_y*cos(0-ref_yaw));
			}
			
			// create a spline
			tk::spline s;
			
			// set (xy) points to the spline
			s.set_points(ptsx, ptsy);
			
			// Define the actual set of points that path planner is going to use
			vector<double> next_x_vals;
          	vector<double> next_y_vals;
            
			// Start with all of the previous path points from last time
			for (int i=0; i < previous_path_x.size(); i++) 
			{
				next_x_vals.push_back(previous_path_x[i]);
				next_y_vals.push_back(previous_path_y[i]);
			}	
			
			// Calculate how to breakup the spline points for the distances(30m) so that we can travel at desired speed
			// We know each distance between points is .2m
			// target velocity 49. mph
			// N * .2m * 49.5mph = (distance between refence and 30m points)
			double target_x = 30.0;
			double target_y = s(target_x);
			double target_dist = sqrt((target_x*target_x)+(target_y*target_y));
			
			double x_add_on = 0;
			
			// Always add 50 points to next_x_vals inclusing previous path points
			for (int i=1; i < 50-previous_path_x.size(); i++) 
			{
			  	double N = (target_dist/(.02*ref_vel/2.24));
				double x_point = x_add_on+(target_x)/N;
				double y_point = s(x_point);
				
				x_add_on = x_point;
				
				double x_ref = x_point;
				double y_ref = y_point;
				
				// Going back from local to Global corordinates
				// rotate back to normal after rotating it earlier
				x_point = (x_ref * cos(ref_yaw)-y_ref*sin(ref_yaw));
				y_point = (x_ref * sin(ref_yaw)+y_ref*cos(ref_yaw));
				
				x_point += ref_x;
				y_point += ref_y;
				
				next_x_vals.push_back(x_point);
				next_y_vals.push_back(y_point);
			}	
			
			

          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          	//double dist_inc = 0.2;
          	// for (int i = 0; i < 50; i++) {
                // double next_s = car_s + (i+1)*dist_inc;
                // double next_d = 6; // 6m
                // vector<double> xy = getXY(next_s, next_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
                // next_x_vals.push_back(xy[0]);
                // next_y_vals.push_back(xy[1]);
          	// }
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
