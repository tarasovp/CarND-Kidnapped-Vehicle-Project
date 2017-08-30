/*
 * map.h
 *
 *  Created on: Dec 12, 2016
 *      Author: mufferm
 */

#ifndef MAP_H_
#define MAP_H_

struct LandmarkObs {
    
    int id;				// Id of matching landmark in the map.
    double x;			// Local (vehicle coordinates) x position of landmark observation [m]
    double y;			// Local (vehicle coordinates) y position of landmark observation [m]
};


class Map {
public:
	
	/*struct single_landmark_s{

		int id_i ; // Landmark ID
		float x_f; // Landmark x-position in the map (global coordinates)
		float y_f; // Landmark y-position in the map (global coordinates)
	};*/

	std::vector<LandmarkObs> landmark_list ; // List of landmarks in the map

};



#endif /* MAP_H_ */
