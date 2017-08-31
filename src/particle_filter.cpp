/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "particle_filter.h"
#include "helper_functions.h"
#define eps 0.0001
using namespace std;

// declare a random engine to be used across multiple and various method calls
static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
    // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
    //   x, y, theta and their uncertainties from GPS) and all weights to 1.
    // Add random Gaussian noise to each particle.
    // NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    
    num_particles = 100;
    weights =  vector<double>(num_particles);
    
    // This line creates a normal (Gaussian) distribution for x
    normal_distribution<double> dist_x(x, std[0]);
    // Create normal distributions for y and theta
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);
    
    for (int i = 0; i < num_particles; ++i) {
        // TODO: Sample  and from these normal distrubtions like this:
        //	 sample_x = dist_x(gen);
        //	 where "gen" is the random engine initialized earlier.
        
        Particle p;
        p.id = i;
        p.x = dist_x(gen);
        p.y = dist_y(gen);
        p.theta = dist_theta(gen);
        p.weight = 1.0;
        
        particles.push_back(p);
        
    }
    
     is_initialized = true;
    
}


void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    // TODO: Add measurements to each particle and add random Gaussian noise.
    // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
    //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    //  http://www.cplusplus.com/reference/random/default_random_engine/
    
    
    // This line creates a normal (Gaussian) distribution for x
    normal_distribution<double> dist_x(0, std_pos[0]);
    // Create normal distributions for y and theta
    normal_distribution<double> dist_y(0, std_pos[1]);
    normal_distribution<double> dist_theta(0, std_pos[2]);
    
    //constatns for all particles
    double yaw_rate_t =yaw_rate*delta_t;
    double v_delta_t = velocity * delta_t;
    double v_yaw_rate;
    if (yaw_rate>eps){v_yaw_rate = velocity / yaw_rate;}
    
    
    for (int i = 0; i < num_particles; ++i) {
        if (yaw_rate>eps)
        {
            
            particles[i].x += v_yaw_rate * (sin(particles[i].theta + yaw_rate_t) - sin(particles[i].theta));
            particles[i].y += v_yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate_t));
        }
        else
        {
            particles[i].x += v_delta_t * cos(particles[i].theta);
            particles[i].y += v_delta_t * sin(particles[i].theta);
        }
        
        particles[i].x+=dist_x(gen);
        particles[i].y+=dist_y(gen);
        particles[i].theta+= yaw_rate_t +dist_theta(gen);
        
    }
    
    
    
    
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
    // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
    //   observed measurement to this particular landmark.
    // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
    //   implement this method and use it as a helper during the updateWeights phase.
    
    
    
    for (int i = 0; i < observations.size(); ++i) {
        double d = numeric_limits<double>::max();
        int id = -1;
        
        //argmin for each measurement
        for (int j = 0; j < predicted.size(); ++j) {
            double now = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
            if (now < d) {
                d = now;
                //number in list - not id!
                id=j;
            }
        }
        
        observations[i].id = id;
        
    }
    
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   std::vector<LandmarkObs> observations, Map map_landmarks) {
    // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
    //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
    //   according to the MAP'S coordinate system. You will need to transform between the two systems.
    //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
    //   The following is a good resource for the theory:
    //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
    //   and the following is a good resource for the actual equation to implement (look at equation
    //   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account
    //   for the fact that the map's y-axis actually points downwards.)
    //   http://planning.cs.uiuc.edu/node99.html
    
    
    //constants for weight transformation
    double a = ( 1/(2*M_PI*std_landmark[0]*std_landmark[1]));
    double b= (2*pow(std_landmark[0], 2));
    double c= (2*pow(std_landmark[1], 2));
    
    for (int i = 0; i < num_particles; ++i) {
        vector<LandmarkObs> landmarks;
        
        // sensor range filter
        for (int j = 0; j <  map_landmarks.landmark_list.size(); ++j)
            if (range_filter(particles[i], map_landmarks.landmark_list[j],sensor_range))
                landmarks.push_back(map_landmarks.landmark_list[j]);
        
        // coordinat transormation
        vector<LandmarkObs> t_observations;
        for (int j = 0; j <  observations.size(); ++j)
            t_observations.push_back (coordinate_transform( particles[i], observations[j]));
        
        
        // dataAssociation
        dataAssociation(landmarks, t_observations);
        
        particles[i].weight = 1.0;
        
        for (unsigned int j = 0; j < t_observations.size(); j++)
            particles[i].weight *= GetWeight (t_observations[j],landmarks[t_observations[j].id],a,b,c);
        
        weights[i]=particles[i].weight;
    }
}

void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional to their weight.
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    
    vector<Particle> new_particles (num_particles);
    discrete_distribution<int> d(weights.begin(), weights.end());
    
    for (int i = 0; i < num_particles; ++i)
        new_particles[i] = particles[d(gen)];
    
    particles = new_particles;
   
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
