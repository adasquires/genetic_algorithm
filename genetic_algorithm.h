#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <pthread.h>
#include <random>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

#include "../consts.h"

#include "../modules/Collide.h"
#include "../modules/FitnessCalc.h"
#include "../modules/TSearch.h"
#include "../modules/util.h"
#include "../modules/VectorMatrix.h"
#include "../modules/Worm.h"

#include "../modules/packages/cxxopts.hpp"
#include "../modules/packages/json.hpp"

#pragma once

#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

using nlohmann::json;

WormBody b; RandomState rs;
double angle = 1.437;
double t_food_start = 1.0;

double angle; double t_food_start = 1.0; double duration = 100.0;

const int gen_count = 300; const int individuals = 60;

double mutation_rate = 0.005; double crossover_rate = 0.7; int tournament_size = individuals / 20; int elitism_size = individuals / 3;

std::array<std::map<std::string, double>, individuals> crossed;

double alpha; double beta; double foodPos_x; double foodPos_y; double p_gamma; double kappa; double lambda;

double AWA_AIY; double AIY_AIY; double AIY_RIA;
double RIA_RMDD; double RIA_RMDV; double SMDD_RIA; double SMDV_RIA; double RIA_RIA;
double RIA_SMDD; double RIA_SMDV; double RMDD_RIA; double RMDV_RIA;
double SMDD_SMDD; double SMDV_SMDV; double RMDD_RMDD; double RMDV_RMDV;
double SMDD_SMDV; double SMDV_SMDD; double SMDD_RMDV; double SMDV_RMDD;
double RMDD_RMDV; double RMDV_RMDD; double SMDD_RMDD_ele; double SMDV_RMDV_ele ; double RMDV_RMDD_ele;

double AIY_tau; double AIY_theta;
double AWA_tau; double AWA_theta;
double RIA_tau; double RIA_theta;
double RMDD_tau; double RMDD_theta;
double RMDV_tau; double RMDV_theta;
double SMDD_tau; double SMDD_theta;
double SMDV_tau; double SMDV_theta;

double NMJ_DB; double NMJ_DD; double NMJ_RMDD; double NMJ_RMDV; double NMJ_SMDD; double NMJ_SMDV;
double NMJ_VBA; double NMJ_VBP; double NMJ_VDA; double NMJ_VDP;

double sr_headgain; double sr_vcgain;

double DB_DB; double VBA_VBA; double VBP_VBP; double DD_DD; double VDA_VDA; double VDP_VDP; double DB_DD;
double VBA_VDA; double VBP_VDP; double DB_VDA; double DB_VDP; double VBA_DD; double VBP_DD; double DD_VDA;
double DD_VDA_ele; double DD_VDP_ele; double VDA_VDP_ele; double VBA_VBP_ele;

double fwd_DB_DB; double fwd_VBP_VBA; double fwd_DD_DD; double fwd_VDP_VDA; double fwd_VBP_DB;

double DB_tau; double DB_theta; double DD_tau; double DD_theta; double VBA_tau; double VBA_theta;
double VBP_tau; double VBP_theta; double VDA_tau; double VDA_theta; double VDP_tau; double VDP_theta;

double DB_tau; double DB_theta; double DD_tau; double DD_theta; double VBA_tau; double VBA_theta;
double VBP_tau; double VBP_theta; double VDA_tau; double VDA_theta; double VDP_tau; double VDP_theta;

class GeneticAlgorithm {

    public:
        GeneticAlgorithm() {};

    // Set random seed for simulation.
    long set_seed(json & simulation_params)
    {
        long seed = static_cast<long>(time(NULL));
        return seed;
    }

    // Run simulation, return fitness data.
    std::tuple<double, double, double, double, double, std::vector<double> > EvaluationFunction(Worm w,
                               			                                        RandomState &rs,
                               			                                        double angle,
                               			                                        std::vector<CollisionObject> & collObjs,
                               			                                        double t_food_start,
                               			                                        std::string output_file)
    {

    std::cout << std::setprecision(10);

    std::ofstream output(output_file);

    std::string output_dir("data/run/NULL/");

    std::ofstream fitfile;
    fitfile.open(output_dir + "fitness.yaml");
    std::ofstream bodyfile;
    bodyfile.open(output_dir + "body.dat");
    std::ofstream actfile;
    actfile.open(output_dir + "act.dat");
    std::ofstream curvfile;
    curvfile.open(output_dir + "curv.dat");

    TVector<double> curvature(1, N_curvs);
    TVector<double> antpostcurv(1, 2);
    antpostcurv.FillContents(0.0);

    w.InitializeState(rs, angle, collObjs);

    FitnessCalc fcalc(w);

    std::vector<double> fitness_concentration;
    double closest = 100.0;
    double closest_time = 0.0;

    double xtp; double ytp; double xt; double yt;

    for (double t = 0.0; t <= DURATION; t += STEPSIZE)
    {
        if (t > t_food_start)
        {
            w.chemo_re.enabled = true;
        }

        w.Step(STEPSIZE, 1);

        fcalc.update();

        w.Curvature(curvature);
        curvfile << curvature << std::endl;

        w.DumpBodyState(bodyfile, skip);

        w.DumpActState(actfile, skip);

        xtp = xt; ytp = yt;
        xt = w.CoMx(); yt = w.CoMy();

        double prox = std::sqrt((xt - foodPos_x) * (xt - foodPos_x) + (yt - foodPos_y) * (yt - foodPos_y));
        if (prox < closest) {
            closest = prox;
            closest_time = t;
        }

        double bodyorientation = w.Orientation();
        double movementorientation = atan2(yt-ytp,xt-xtp);
        double anglediff = movementorientation - bodyorientation;
        double concentration = w.chemo_re.get_concentration(VecXY(w.b.X(1), w.b.Y(1)));
        fitness_concentration.push_back(concentration);

    }

    bodyfile.close();
    actfile.close();
    curvfile.close();
    fitfile << fcalc.strprintf();

    double distance_travelled = fcalc.distancetravelled;

    std::cout << "fcalc distance travelled: " << distance_travelled << std::endl;

    return std::make_tuple(distance_travelled, w.CoMx(), w.CoMy(), closest, closest_time, fitness_concentration);

    }

    // Return random double for parameter randomization.
    double get_rand(double min,
                    double max)
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(min, max);

        double rand_num = dis(gen);

        while (rand_num == 0.0)
        {
            rand_num = dis(gen);
        }

        return rand_num;
    }

    // Return fitness.
    std::tuple<double, double, double> fitness(double distance_travelled,
                                               double head_x,
                                               double head_y,
                                               double closest,
                                               double closest_time,
                                               std::vector<double> fitness_concentration,
                                               std::string input)
    {

        std::cout << "fit distance travelled: " << distance_travelled << std::endl;

        std::ifstream f(input);

        if (!f)
        {
            std::cerr << "params.json not found" << std::endl;
   	    }

        json data;
        f >> data;

        // Proximity fitness.

        foodPos_x = data["ChemoReceptors"]["foodPos"]["x"].get<double>();
        foodPos_y = data["ChemoReceptors"]["foodPos"]["y"].get<double>();

        double dx = head_x - foodPos_x;
        double dy = head_y - foodPos_y;

        double max_dist = std::sqrt(foodPos_x * foodPos_x + foodPos_y * foodPos_y);
        //double max_dist = AvgSpeed * DURATION;
        // << "max dist: " << max_dist << std::endl;
        double dist_to_food = std::sqrt(dx*dx + dy*dy);
        //std::cout << "dist to food: " << dist_to_food << std::endl;
        double close = std::abs((closest / max_dist) * (closest_time / (DURATION * STEPSIZE)));
        //std::cout << "close: " << close << std::endl;
        double raw_proximity = 1.0 - dist_to_food / max_dist;
        double proximity = raw_proximity * (1.0 - close);
        proximity = std::max(0.0, std::min(1.0, proximity));
        //std::cout << "Proximity: " << proximity << std::endl;

        // Velocity fitness.

        double v = distance_travelled / DURATION;
        std::cout << "Distance travelled: " << distance_travelled << std::endl;
        //double velocity = 1 - (std::abs(AvgSpeed - v) / AvgSpeed);
        // velocity = std::max(0.0, std::min(1.0, velocity));
        double diff = std::abs((v - AvgSpeed) / AvgSpeed + 1e-6);
        double velocity = 1.0 / (1.0 + diff);
        if (!std::isfinite(velocity)) velocity = 0.0;
        std::cout << "Velocity: " << velocity << std::endl;

        // Chemosensation fitness.

        double chemosensation_fitness = 0;
        size_t count = 0;
        for (size_t i = 1; i < fitness_concentration.size(); ++i)
        {
            if (fitness_concentration[i] > fitness_concentration[i-1])
            {
                count += 1;
            }
        }

        double chemosensation_count = 0;
        double chemosensation_value = 0;
        for (int i = 0; i < fitness_concentration.size(); ++i) {
            if (fitness_concentration[i] > chemosensation_value) {
                chemosensation_value = fitness_concentration[i];
                chemosensation_count = static_cast<double>(i);
            }
        }

        chemosensation_count /= (duration / STEPSIZE);
        chemosensation_value /= 244.00503969854154;

        //std::cout << "Chemosensation value: " << chemosensation_value << std::endl;
        //std::cout << "Chemosensation count: " << chemosensation_count << std::endl;

        chemosensation_fitness = static_cast<double>(count) / (fitness_concentration.size() - 1);
        chemosensation_fitness = std::max(0.0, std::min(1.0, chemosensation_fitness));

        double chemosensation = 0.5 * chemosensation_count + 0.5 * chemosensation_value;
        //std::cout << "Chemosensation fitness: " << chemosensation << std::endl;

        double cumulative_concentration = std::accumulate(fitness_concentration.begin(), fitness_concentration.end(), 0.0);
        double height = 244.00503969854154;
        double width = (DURATION / STEPSIZE);
        cumulative_concentration /= (0.5 * (height * 0.39 * width) + (height * 0.61 * width));
        cumulative_concentration = std::max(0.0, std::min(1.0, cumulative_concentration));
        std::cout << "Cumulative concentration fitness: " << cumulative_concentration << std::endl;

        std::vector<double> chemosensation_derivative;
        for (int i = 1; i < fitness_concentration.size(); ++i) {
            chemosensation_derivative.push_back(fitness_concentration[i] - fitness_concentration[i-1]);
        }

        for (size_t i = 1; i < fitness_concentration.size(); ++i) {

        }
        double dt = (DURATION / STEPSIZE) / fitness_concentration.size();
        double r0 = max_dist;
        double integral = 0.0;

        for (int i = 0; i <= fitness_concentration.size(); ++i) {
            double t = i * dt;
            double weight = (i == 0 || i == fitness_concentration.size()) ? 0.5 : 1.0;
            integral += weight * (fitness_concentration[i] / r0);
        }

        integral *= dt;

        double chemotaxis_index = 1.0 - (1.0 / (DURATION / STEPSIZE)) * integral;

        //std::cout << "Chemotaxis: " << chemotaxis_index << std::endl;

        double chemo_dist = exp(fitness_concentration.back() * fitness_concentration.back()) / (2.0 * (max_dist * max_dist));
        //std::cout << "chemo dist: " << chemo_dist << std::endl;

        double chemo_fit = fitness_concentration.back() / 244.00503969854154;

        return std::make_tuple(velocity, proximity, cumulative_concentration);

    };

    // Return indices of fittest individuals.
    std::vector<int> get_fittest(std::vector<double> arr)
	{
    	// Ignore NaN values.
    	std::vector<int> valid_indices;
    	for (int i = 0; i < arr.size(); ++i)
    	{
        	if (!std::isnan(arr[i]))
        	{
            	valid_indices.push_back(i);
        	}
    	}

    	if (valid_indices.size() <= individuals)
    	{
        	std::sort(valid_indices.begin(),
                  	valid_indices.end(),
                  	[&](int i, int j)
                    { return arr[i] > arr[j]; });
        	return valid_indices;
    	}

    	// Partial sort to get the top individuals.
    	std::partial_sort(valid_indices.begin(),
                      	valid_indices.begin() + individuals,
                      	valid_indices.end(),
                      	[&](int i, int j)
                        { return arr[i] > arr[j]; });

   		valid_indices.resize(individuals);

    	return valid_indices;
	}

    // Hamming distance (measure of generation diversity).
    double hamming_distance(std::array<std::map<std::string, double>, individuals>& params)
    {
        int distance = 0;
        size_t num_params = params[0].size();

        for (size_t i = 0; i < individuals; ++i)
        {
            for (size_t j = i + 1; j < individuals; ++j)
            {
                for (const auto& [key, value1] : params[i])
                {
                    double value2 = params[j].at(key);
                    if (std::abs(value1 - value2) > 1e-6)
                    {
                        ++distance;
                    }
                }
            }
        }

        size_t num_pairs = individuals * (individuals - 1) / 2;
        size_t max_possible_differences = num_pairs * num_params;

        double d = static_cast<double>(distance) / static_cast<double>(max_possible_differences);
        d = std::max(0.0, std::min(1.0, d));

        return d;
    }


    std::map<std::string, double> elitism_select(std::array<std::map<std::string, double>, individuals>& params,
                                                 std::vector<double>& fitness_values,
                                                 int n)
    {
        std::vector<int> fitness = get_fittest(fitness_values);
        return params[fitness[n]];
    }

    // Fitness proportional selection.
    std::map<std::string, double> fitness_proportional_select(
            const std::array<std::map<std::string, double>, individuals>& gen_params,
            const std::vector<double>& fitness)
    {
        double epsilon = 1e-6; // avoid zero probabilities
        double total_fitness = 0.0;
        std::vector<double> adjusted_fitness(fitness.size());

        for (size_t i = 0; i < fitness.size(); i++) {
            adjusted_fitness[i] = fitness[i] + epsilon;
            total_fitness += adjusted_fitness[i];
        }

        // Create cumulative distribution
        std::vector<double> cumulative(fitness.size(), 0.0);
        cumulative[0] = adjusted_fitness[0] / total_fitness;
        for (size_t i = 1; i < fitness.size(); i++) {
            cumulative[i] = cumulative[i-1] + adjusted_fitness[i] / total_fitness;
        }

        // Draw random number
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        double r = dist(gen);

        // Find selected index
        size_t selected_index = 0;
        for (; selected_index < cumulative.size(); selected_index++) {
            if (r <= cumulative[selected_index]) break;
        }

        return gen_params[selected_index];
    }

    // Linear rank selection.
    std::map<std::string, double> linear_rank_select(std::array<std::map<std::string, double>, individuals>& params,
                                                     std::vector<double>& fitness_values,
                                                     double n)
    {

        double n_plus = 2 - n;
        double n_minus = n;

        // Rank individuals from worst to best based on fitness.

        std::array<std::map<std::string, double>, individuals / 2> selected;
        std::vector<int> fitness = get_fittest(fitness_values);

        std::vector<int> ranked = fitness;
        std::reverse(ranked.begin(), ranked.end());

        // Get probability based on linear rank function.

        std::vector<double> probabilities;
        probabilities.reserve(fitness.size());
        for (size_t i = 1; i < ranked.size() + 1; ++i)
        {
            std::cout << i << std::endl;
            double probability = ((1.0/individuals) * (n_minus + (n_plus - n_minus) * ((i - 1.0) / (individuals - 1.0))));
            probabilities.push_back(probability);
        }

        // Roulette wheel.

        double fitness_sum = std::accumulate(probabilities.begin(), probabilities.end(), 0.0);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0, fitness_sum);

        double cumulative_sum = 0.0;
        double random_number = dis(gen);
        for (size_t i = 0; i < ranked.size(); ++i)
        {
            cumulative_sum += probabilities[i];
            if (random_number <= cumulative_sum)
            {
                return params[ranked[i]];
        	}
		}
    }

    // Split rank selection.
    std::array<std::map<std::string, double>, individuals / 2> split_rank_select(std::array<std::map<std::string, double>, individuals>& params,
                                                                                 std::vector<double>& fitness_values,
                                                                                 double select_lambda)
    {

        // Rank indivdiuals from worst to best based on fitness. .

        std::array<std::map<std::string, double>, individuals / 2> selected;
        std::vector<int> fitness = get_fittest(fitness_values);

        std::vector<int> ranked = fitness;
        std::reverse(ranked.begin(), ranked.end());

        // Scale to probability based on split rank function.

        std::vector<double> probabilities;
        for (size_t i = 1; i < ranked.size() + 1; ++i)
        {
            if (i <= individuals / 2)
            {
            	double probability = (1 - select_lambda) * ( (8.0 * i) / (individuals * (individuals + 2.0)) );
            	probabilities.push_back(probability);
            } else
            {
              	double probability = (1 - select_lambda) * ( (8.0 * i) / (3.0*individuals * (individuals + 2.0)) );
            	probabilities.push_back(probability);
            }

        }

        // Roulette wheel.

        double fitness_sum = std::accumulate(probabilities.begin(), probabilities.end(), 0.0);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0, fitness_sum);

        double cumulative_sum = 0.0;
        int j = 0;
        while (j < individuals / 2)
        {
            double random_number = dis(gen);
        	for (size_t i = 0; i < ranked.size(); ++i)
        	{
                cumulative_sum += probabilities[i];
            	if (random_number < cumulative_sum)
            	{
                	selected[j] = params[ranked[i]];
                	j++;
                    break;
            	}
        	}
		}
        return selected;
    }

    // Tournament selection.
    std::array<std::map<std::string, double>, individuals / 2> tournament_select(std::array<std::map<std::string, double>, individuals>& params,
                                                                                 std::vector<double>& scaled_fitness,
                                                                                 int tournament_size,
                                                                                 int elitism_size)
    {

        // Get fitness values.
        std::array<std::map<std::string, double>, individuals / 2> selected;
        std::vector<int> fitness = get_fittest(scaled_fitness);

        // Elitism to always select the fittest individuals.
        for (int i = 0; i < elitism_size; ++i)
        {
            selected[i] = params[fitness[i]];
        }

        // Tournament selection.

        std::vector<int> used;

        for (int i = elitism_size; i < (individuals / 2); ++i)
        {
    	    int k = rand() % (individuals / 2);

    	    for (int i = 0; i < tournament_size; ++i)
    	    {
        	    int j = rand() % (individuals / 2);
                auto it = std::find(used.begin(), used.end(), j);
        	    if (fitness[j] > fitness[k] && it == used.end())
        	    {
            	    k = j;
        	    }
    	    }

            selected[i] = params[fitness[k]];
            used.push_back(k);

        }

        return selected;
    }

    // Adaptive crossover to generate offspring.
    std::map<std::string, double> crossover(std::map<std::string, double> params_1,
                                            std::map<std::string, double> params_2,
                                                std::vector<std::string> param_names,
                                                double crossover_rate)
    {

        std::map<std::string, double> child;

        if ((rand() % 100) / 100.0 < crossover_rate)
        {
            int crossover_point = std::rand() % param_names.size();
            for (int i = 0; i < crossover_point; ++i)
            {
                std::string p = param_names[i];
                std::cout << p << std::endl;
                child[p] = params_1.at(p);
            }
            for (int i = crossover_point; i < params_1.size(); ++i)
            {
                std::string p = param_names[i];
                std::cout << p << std::endl;
                child[p] = params_2.at(p);
            }
        } else
        {
            if ((std::rand() % 2) + 1 == 1)
            {
                for (int i = 0; i < param_names.size(); ++i)
                {
                    std::string p = param_names[i];
                    std::cout << p << std::endl;
                    child[p] = params_1.at(p);
                }
            } else
            {
                for (int i = 0; i < param_names.size(); ++i)
                {
                    std::string p = param_names[i];
                    std::cout << p << std::endl;
                    child[p] = params_2.at(p);
                }
            }
        }

        return child;

    };

    // Mutate individuals based on Gaussian distribution.
    std::map<std::string, double> mutate(std::map<std::string, double>& params,
                                        double mutation_rate,
                                        double mutation_sigma,
                                        std::vector<std::string> param_names,
                                        std::vector<std::string> change_params,
                                        std::string output_file)
    {
        std::map<std::string, double> mutated_geno;
        mutated_geno = params;
        std::random_device rd;
        std::mt19937 gen(rd());
        //std::normal_distribution<> dist(0.0, mutation_sigma);
        std::uniform_real_distribution<> prob_dist(0.0, 1.0);

        //for (int i = 0; i < change_params.size(); ++i)
        //{
        //    if (get_rand(0.0, 1.0) < mutation_rate)
       	//    {
        //        std::string p = change_params[i];
        //        if (params[p] > 0.0)
        //        {
        //            mutated_geno[p] += std::abs(dist(gen));
        //        } else if (params[p] < 0.0)
        //        {
        //            mutated_geno[p] -= std::abs(dist(gen));
        //        }
        //    }
        //}

        std::uniform_real_distribution<> dist(0.0, 1.0);
        std::normal_distribution<> perturb(0.0, mutation_sigma);

        for (auto& p: params) {
            if (dist(gen) < mutation_rate) {
                p.second += perturb(gen);
            }
        }

        // Ensure dorsal/ventral symmetry.
   	    mutated_geno["RIA_RMDD"] = mutated_geno["RIA_RMDV"];
        mutated_geno["SMDD_RIA"] = mutated_geno["SMDV_RIA"];
        mutated_geno["RIA_SMDD"] = mutated_geno["RIA_SMDV"];
        mutated_geno["RMDD_RIA"] = mutated_geno["RMDV_RIA"];
        mutated_geno["SMDD_SMDD"] = mutated_geno["SMDV_SMDV"];
        mutated_geno["RMDD_RMDD"] = mutated_geno["RMDV_RMDV"];
        mutated_geno["SMDD_SMDV"] = mutated_geno["SMDV_SMDD"];
        mutated_geno["SMDD_RMDV"] = mutated_geno["SMDV_RMDD"];
        mutated_geno["RMDD_RMDV"] = mutated_geno["RMDV_RMDD"];
        mutated_geno["SMDD_RMDD_ele"] = mutated_geno["SMDV_RMDV_ele"];
        mutated_geno["RMDD_theta"] = mutated_geno["RMDV_theta"];
        mutated_geno["SMDD_tau"] = mutated_geno["SMDV_tau"];
        mutated_geno["SMDD_theta"] = mutated_geno["SMDV_theta"];
        mutated_geno["SMDD_RMDD_ele"] = std::abs(mutated_geno["SMDD_RMDD_ele"]);
   	    mutated_geno["SMDV_RMDV_ele"] = std::abs(mutated_geno["SMDV_RMDV_ele"]);
        mutated_geno["RMDV_RMDD_ele"] = std::abs(mutated_geno["RMDV_RMDD_ele"]);
        mutated_geno["foodPos_x"] = 0.002;
        mutated_geno["foodPos_y"] = 0.003;

        return mutated_geno;

    }

};


#endif //GENETIC_ALGORITHM_H
