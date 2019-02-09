/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

namespace {
    double multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs,
                       double mu_x, double mu_y) {
        // calculate normalization term
        double gauss_norm;
        gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

        // calculate exponent
        double exponent;
        exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
                   + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));

        // calculate weight using normalization terms and exponent
        double weight;
        weight = gauss_norm * exp(-exponent);

        return weight;
    }

//    template<class T>
//    T constrainAngle(T x){
//        T pi_2 = 2 * M_PI;
//        x = fmod(x + M_PI,pi_2);
//        if (x < 0)
//            x += pi_2;
//        return x - M_PI;
//    }
}

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 1000;  // TODO: Set the number of particles

  std::default_random_engine gen;
  double std_x, std_y, std_theta;  // Standard deviations for x, y, and theta
  std_x = std[0];
  std_y = std[1];
  std_theta = std[2];

  // This line creates a normal (Gaussian) distribution for x
  std::normal_distribution<double> dist_x(x, std_x);
  std::normal_distribution<double> dist_y(y, std_y);
  std::normal_distribution<double> dist_theta(theta, std_theta);

  for (int i = 0; i < num_particles; ++i) {
    Particle p;
    p.id = i;
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);//constrainAngle(dist_theta(gen));
    p.weight = 1;
    particles.push_back(p);
  }

  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  std::default_random_engine gen;
  double std_x, std_y, std_theta;  // Standard deviations for x, y, and theta
  std_x = std_pos[0];
  std_y = std_pos[1];
  std_theta = std_pos[2];

  for (int i = 0; i < num_particles; ++i) {

      double pred_x, pred_y, pred_theta;
      if(yaw_rate < 1e-10) {
          pred_x = particles[i].x + velocity * delta_t * cos(particles[i].theta);
          pred_y = particles[i].y + velocity * delta_t * sin(particles[i].theta);
          pred_theta = particles[i].theta;
      } else {
          pred_x = particles[i].x + velocity / yaw_rate * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
          pred_y = particles[i].y + velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
          pred_theta = particles[i].theta + yaw_rate * delta_t;
      }
    std::normal_distribution<double> dist_x(pred_x, std_x);
    std::normal_distribution<double> dist_y(pred_y, std_y);
    std::normal_distribution<double> dist_theta(pred_theta, std_theta);

    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);//constrainAngle(dist_theta(gen));
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
   for(int i = 0; i < observations.size(); ++i) {
       auto obs_x = observations[i].x;
       auto obs_y = observations[i].y;

       auto smallest_distance = std::numeric_limits<double>::infinity();
       for(auto& pred : predicted) {
           auto position_rmse = dist(obs_x, obs_y, pred.x, pred.y);
           if(position_rmse < smallest_distance) {
               smallest_distance = position_rmse;
               observations[i].id = pred.id;
           }
       }
   }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
   weights.clear();
   for (int i = 0; i < num_particles; ++i) {
    double x_part = particles[i].x;
    double y_part = particles[i].y;
    double cos_theta = std::cos(particles[i].theta);
    double sin_theta = std::sin(particles[i].theta);

    // computes the inverse of a transformation matrix
    double det = cos_theta * cos_theta +
                 sin_theta * sin_theta;

    double invdet = 1 / det;

    double minv00 = cos_theta * invdet;
    double minv01 =  sin_theta * invdet;
    double minv02 = (-sin_theta * y_part - x_part * cos_theta) * invdet;
    double minv10 =  -sin_theta * invdet;
    double minv11 = cos_theta * invdet;
    double minv12 = (sin_theta * x_part - cos_theta * y_part) * invdet;

    std::vector<LandmarkObs> mapInVeh;
    for (auto& map_landmark : map_landmarks.landmark_list) {
      double distance = dist(x_part, y_part, map_landmark.x_f, map_landmark.y_f);
      if(distance <= sensor_range) {
        // transform to VEHICLE'S coordinate system
        double x_map_in_veh = minv02 + (minv00 * map_landmark.x_f) + (minv01 * map_landmark.y_f);
        double y_map_in_veh = minv12 + (minv10 * map_landmark.x_f) + (minv11 * map_landmark.y_f);
        LandmarkObs landmark;
        landmark.id = map_landmark.id_i;
        landmark.x = x_map_in_veh;
        landmark.y = y_map_in_veh;
        mapInVeh.push_back(landmark);
      }
    }

    // associate obs to landmark
    std::vector<LandmarkObs> observations_ = observations;
    dataAssociation(mapInVeh, observations_);

    // calculate weight
    double total_weight = 1.0;
    std::vector<int> associations;
    std::vector<double> sense_x;
    std::vector<double> sense_y;
    for(auto& obs : observations_) {
        auto it = std::find_if(mapInVeh.begin(), mapInVeh.end(), [obs](LandmarkObs const& obj){ return obj.id == obs.id; });
        if(it != mapInVeh.end()) {
            total_weight *= multiv_prob(std_landmark[0], std_landmark[1], obs.x, obs.y, it->x, it->y);

            auto it_map = std::find_if(map_landmarks.landmark_list.begin(), map_landmarks.landmark_list.end(), [it](Map::single_landmark_s const& obj){ return obj.id_i == it->id; });

            associations.push_back(it_map->id_i);
            sense_x.push_back(it_map->x_f);
            sense_y.push_back(it_map->y_f);
        }
    }

    particles[i].weight = (total_weight < 1.0) ? total_weight : particles[i].weight;

    particles[i].associations.clear();
    particles[i].sense_x.clear();
    particles[i].sense_y.clear();
    SetAssociations(particles[i],associations,sense_x,sense_y);

    // for std::discrete_distribution
    weights.push_back(particles[i].weight);
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
   std::random_device rd;
   std::mt19937 gen(rd());
   std::discrete_distribution<> d(weights.begin(),weights.end());
   //double total_weight = std::accumulate(weights.begin(),weights.end(),0.0);
   std::vector<Particle> resampled_particles;
   for (int i = 0; i < num_particles; ++i) {
       Particle resampled_part = particles[d(gen)];
       resampled_particles.push_back(resampled_part);
   }
   particles.clear();
   particles = resampled_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}