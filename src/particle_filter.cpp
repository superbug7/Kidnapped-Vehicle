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
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	if(!is_initialized)
	{
		num_particles = 15;
		
		default_random_engine gen;

		normal_distribution<double> dist_x(x, std[0]);
  		normal_distribution<double> dist_y(y, std[1]);
  		normal_distribution<double> dist_theta(theta, std[2]);

		for(int i =0; i< num_particles;i++)
		{
			printf("initializing %d \n", i);
			Particle P;
			P.id= i;
			P.x = dist_x(gen);
			P.y = dist_y(gen);
			P.theta = dist_theta(gen);
			P.weight= 1.0;
			particles.push_back(P);
		}
	}
	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;
	double x_new;
	double y_new;
	double theta_new;

	for( int i=0 ; i< num_particles; i++)
	{

		if ( fabs(yaw_rate) < 0.0001){
		x_new = particles[i].x + (velocity * delta_t * cos(theta_new));
		y_new = particles[i].y + (velocity * delta_t * sin(theta_new));
		theta_new = particles[i].theta;} 
		else{
		x_new = particles[i].x + ((velocity/yaw_rate) * (sin(theta_new) - sin(particles[i].theta)));
		y_new = particles[i].y + ((velocity/yaw_rate) * (cos(particles[i].theta - cos(theta_new) )));
		theta_new = particles[i].theta + yaw_rate * delta_t;}
		

		normal_distribution<double> dist_x(x_new, std_pos[0]);
    		normal_distribution<double> dist_y(y_new, std_pos[1]);
    		normal_distribution<double> dist_theta(theta_new, std_pos[2]);


		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);

	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	 

         for(int i =0; i < num_particles;i++)
	 {

		vector<LandmarkObs> obs_transformed;
		vector<LandmarkObs> landmark_predicted;
		double p_x = particles[i].x;
		double p_y = particles[i].y;
		double p_theta = particles[i].theta;

		for(int j =0; j < map_landmarks.landmark_list.size() ;j++)
		{
			Map::single_landmark_s landmark_curr = map_landmarks.landmark_list[j];
			double new_dist = dist(p_x, p_y, landmark_curr.x_f, landmark_curr.y_f);
			if ( new_dist < sensor_range )
				landmark_predicted.push_back( LandmarkObs { landmark_curr.id_i, landmark_curr.x_f, landmark_curr.y_f});
		}

		for( int j =0 ; j < observations.size(); j++)
		{
			LandmarkObs obs_trans;
			obs_trans.id = j;
			obs_trans.x = p_x + (cos(p_theta) * observations[j].x) - (sin(p_theta) *observations[j].y);
			obs_trans.y = p_y + (sin(p_theta) * observations[j].x) + (cos(p_theta) *observations[j].y);	
			obs_transformed.push_back(obs_trans);
		}

		for ( int i=0; i < obs_transformed.size() ; i++)
		{
			double min_dist;
			for ( int j =0; j < landmark_predicted.size(); j++)
			{
				double dista = dist(obs_transformed[i].x, obs_transformed[i].y, landmark_predicted[j].x , landmark_predicted[j].y);
			        if( j == 0)
				{
					min_dist = dista;
					obs_transformed[i].id = landmark_predicted[i].id;
				}
				else if ( dista < min_dist )
				{
					min_dist = dista;
					obs_transformed[i].id = landmark_predicted[i].id;
					 
				}
			}
		}
	
		particles[i].weight =1.0;

		for( int m =0 ; m < obs_transformed.size(); m++)
		{
			for ( int j=0; j < landmark_predicted.size();j++)
			{
				double x = pow(obs_transformed[m].x - landmark_predicted[j].x, 2) / ( 2 * pow(std_landmark[0],2));
				double y = pow(obs_transformed[m].y - landmark_predicted[j].y, 2) / ( 2 * pow(std_landmark[1],2));
				double w = exp(-(x + y )) / ( 2 * M_PI * std_landmark[0] * std_landmark[1]);
				particles[i].weight *= w;
			}
		}				

		weights.push_back(particles[i].weight);		
	}	
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	std::random_device rand;
	std::mt19937 gen(rand());
	std::discrete_distribution<int> weight_distribution(weights.begin(), weights.end());

	vector<Particle> new_particles;

	for( int i =0; i < num_particles ; i++)
	{	
		int idx =weight_distribution(gen);
		new_particles.push_back(particles[idx]);
	}
	particles = new_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

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
