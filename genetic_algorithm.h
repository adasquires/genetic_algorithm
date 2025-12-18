#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <math.h>
#include "nlopt.hpp"
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

double angle; double t_food_start = 1.0; double duration = 100.0;

const int gen_count = 100; const int individuals = 60;

double mutation_rate = 0.005; double crossover_rate = 0.7; int tournament_size = individuals / 10; int elitism_size = individuals / 10;

double random_double;

std::array<json, individuals> gen_params;
std::array<json, individuals/2> selected_params;
std::array<json, individuals> crossed;
std::array<json, individuals> mutated;

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
    std::tuple<double, double> EvaluationFunction(Worm w,
                                                  RandomState &rs,
                                                  double angle,
                                                  std::vector<CollisionObject> & collObjs,
                                                  double t_food_start,
                                                  std::string output_file,
                                                  double foodpos_x,
                                                  double foodpos_y)
    {

        ofstream fitfile;
        fitfile.open("fitnes.yml");

        ofstream bodyfile;
        bodyfile.open("body.dat");

        ofstream actfile;
        actfile.open("act.dat");
        w.DumpActState_header(actfile);

        ofstream curvfile;
        curvfile.open("curv.dat");

        TVector<double> curvature(1, N_curvs);
        TVector<double> antpostcurv(1, 2);
        antpostcurv.FillContents(0.0);

        std::cout << std::setprecision(10);

        w.InitializeState(rs, angle, collObjs);

        FitnessCalc fcalc(w);

        std::vector<double> fitness_concentration;
        double xtp; double ytp; double xt; double yt;
        double distance = 0.0;
        double t,fitness=0.0;
        double fA = 0.0;
        double accdist = 0.0;
        double totaldist = 0.0;
        double MaxDist = dist(VecXY(foodpos_x, foodpos_y), VecXY(0.0, 0.0));
        std::cout << "Food placed at " << foodpos_x << ", " << foodpos_y << std::endl;
        std::cout << "MaxDist: " << MaxDist << std::endl;

        // Time loop
        for (double t = 0.0; t <= DURATION; t += STEPSIZE)
        {
            if (t > t_food_start)
            {
                w.chemo_re.enabled = true;
            }
            // do the actual step
            w.Step(STEPSIZE, 1);
            // update fitness
            fcalc.update();
            w.Curvature(curvature);
            double concentration = dist(VecXY(w.b.X(1), w.b.Y(1)), VecXY(foodpos_x, foodpos_y));
            accdist += concentration;
            fitness_concentration.push_back(concentration);
            distance += dist_sqrd(VecXY(w.b.X(1), w.b.Y(1)), VecXY(xt, yt));
            xt = w.b.X(1);
            yt = w.b.Y(1);
        }

        std::cout << "accdist: " << accdist << std::endl;
        totaldist = (accdist/(DURATION/STEPSIZE));
        std::cout << "totaldist: " << totaldist << std::endl;
        fA = (MaxDist - totaldist)/MaxDist;
        //fA = totaldist / (MaxDist/AvgSpeed) / 2;
        fA = fA < 0 ? 0.0 : fA;
        fA = fA > 1 ? 0.0 : fA;
        std::cout << "fA: " << fA << std::endl;
        fitness += fA;
        std::cout << "Chemotaxis index: " << (MaxDist - totaldist)/MaxDist << std::endl;

        //std::cout << "fcalc distance travelled: " << dist << std::endl;

        w.Curvature(curvature);
        curvfile << curvature << endl;
        w.DumpBodyState(bodyfile, skip);
        w.DumpActState(actfile, skip);

        bodyfile.close();

        actfile.close();

        curvfile.close();

        fitfile << fcalc.strprintf();


        return std::make_tuple(distance, fitness);

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
    std::tuple<double, double> fitness(double distancetravelled,
                                       double fitness_concentration)
    {

        std::cout << "Distance travelled: " << distancetravelled << std::endl;
        if (std::isnan(distancetravelled))
        {
            return std::make_tuple(0.0, 0.0);
        }
        // Velocity fitness.
        double avgdist = AvgSpeed * DURATION;
        double velocity = 1 - (std::abs(avgdist - distancetravelled) / avgdist);
        velocity = std::max(0.0, std::min(1.0, velocity));
        std::cout << "Velocity: " << velocity << std::endl;

        // Chemosensation fitness.
        //double h0 = fitness_concentration[0];
        // compute integral of h(t)/h(0)
        //double integral = 0.0;
        //for (size_t i = 1; i < fitness_concentration.size(); ++i)
        //{
        //    double y1 = fitness_concentration[i-1] / h0;
        //    double y2 = fitness_concentration[i]   / h0;
        //    integral += (0.5 * (y2 + y1)) * STEPSIZE;   // trapezoidal area
        //}

        //double chemotaxis_index = 1.0 - integral / DURATION;
        //chemotaxis_index = std::max(0.0, std::min(1.0, chemotaxis_index));
        //std::cout << "Chemotaxis index: " << chemotaxis_index << std::endl;

        return std::make_tuple(velocity, fitness_concentration);

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
    double hamming_distance(std::array<json, individuals>& params)
    {
        std::vector<json> GenomeI, GenomeJ;
        //std::vector<int> change = {11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 108, 109, 110, 111, 112, 113};
        std::vector<int> change = {0, 1, 4, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63, 67, 71, 75, 79, 83, 87, 91, 95, 99, 103};
        int distance = 0;
        int size = 0;
        for (size_t i = 0; i < individuals; i++)
        {
            GenomeI.clear();
            flatten(params[i], GenomeI);
            for (size_t j = i + 1; j < individuals; j++)
            {
                GenomeJ.clear();
                flatten(params[j], GenomeJ);
                for (int k = 0; k < change.size(); k++)
                {
                    size++;
                    if (std::abs(GenomeI[change[k]].get<double>() - GenomeJ[change[k]].get<double>()) > 1e-3) distance++;
                }
            }
        }

        double d = static_cast<double>(distance) / size;
        d = std::max(0.0, std::min(1.0, d));

        return d;
    }


    json elitism_select(std::array<json, individuals>& params,
                        std::vector<double>& fitness_values,
                        int n)
    {
        std::vector<int> fitness = get_fittest(fitness_values);
        return params[fitness[n]];
    }

    // Fitness proportional selection.
    json fitness_proportional_select(std::array<json, individuals>& params,
                                     std::vector<double>& fitness_values)
    {

        // Sort individuals by fitness.
        std::array<json, individuals / 2> selected;
        std::vector<int> fitness = get_fittest(fitness_values);

        // Roulette wheel.
        double fitness_sum = std::accumulate(fitness_values.begin(), fitness_values.end(), 0.0);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, fitness_sum);

        double rand_val = dis(gen);
        double cum_sum = 0.0;

        for (size_t i = 0; i < fitness_values.size(); ++i)
        {
            cum_sum += fitness_values[i];
            if (rand_val <= cum_sum)
            {
                return params[i];
            }
        }
        return params[fitness[0]];

    }

    // Linear rank selection.
    json linear_rank_select(std::array<json, individuals>& params,
                            std::vector<double>& fitness_values,
                            double n)
    {

        double n_plus = 2 - n;
        double n_minus = n;

        // Rank individuals from worst to best based on fitness.

        std::array<json, individuals / 2> selected;
        std::vector<int> fitness = get_fittest(fitness_values);

        std::vector<int> ranked = fitness;
        std::reverse(ranked.begin(), ranked.end());

        // Get probability based on linear rank function.

        std::vector<double> probabilities;
        for (size_t i = 1; i < ranked.size() + 1; ++i)
        {
            //std::cout << i << std::endl;
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
            if (random_number < cumulative_sum)
            {
                return params[ranked[i]];
            }
        }
    }

    // Split rank selection.
    std::array<json, individuals / 2> split_rank_select(std::array<json, individuals>& params,
                                                        std::vector<double>& fitness_values,
                                                        double select_lambda)
    {

        // Rank indivdiuals from worst to best based on fitness. .

        std::array<json, individuals / 2> selected;
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
    std::array<json, individuals / 2> tournament_select(std::array<json, individuals>& params,
                                                        std::vector<double>& scaled_fitness,
                                                        int tournament_size,
                                                        int elitism_size)
    {

        // Get fitness values.
        std::array<json, individuals / 2> selected;
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

    void flatten(const json& j, std::vector<json>& out) {
        if (j.is_primitive()) {
            out.push_back(j);
        }
        else if (j.is_array()) {
            for (const auto& v : j)
                flatten(v, out);
        }
        else if (j.is_object()) {
            for (const auto& [k, v] : j.items())
                flatten(v, out);
        }
    }

    void unflatten(const json& structure, const std::vector<json>& genome, int& idx, json& out) {
        if (structure.is_primitive()) {
            out = genome[idx++];
        }
        else if (structure.is_array()) {
            out = json::array();
            for (const auto& v : structure) {
                json child;
                unflatten(v, genome, idx, child);
                out.push_back(child);
            }
        }
        else if (structure.is_object()) {
            out = json::object();
            for (const auto& [k, v] : structure.items()) {
                json child;
                unflatten(v, genome, idx, child);
                out[k] = child;
            }
        }
    }

    // Adaptive crossover to generate offspring.
    json crossover(const json& p1, const json& p2, double rate) {
        // flatten parents
        std::vector<json> g1, g2;
        flatten(p1, g1);
        flatten(p2, g2);

        if (g1.size() != g2.size()) {
            throw std::runtime_error("Parents have different genome lengths.");
        }

        int N = g1.size();
        std::vector<json> child_genome(N);

        // no crossover? random parent
        if (((rand() % 100) / 100.0) >= rate) {
            return (rand() % 2 == 0 ? p1 : p2);
        }

        // choose 1-point
        int pt = rand() % N;

        for (int i = 0; i < N; i++) {
            child_genome[i] = (i < pt ? g1[i] : g2[i]);
        }

        // rebuild JSON from genome
        json child;
        int idx = 0;
        unflatten(p1, child_genome, idx, child);

        return child;
    }


    // Mutate individuals based on Gaussian distribution.
    json mutate(json& params,
                double mutation_rate,
                double mutation_sigma)
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<> dist(0.0, mutation_sigma);
        std::uniform_real_distribution<> prob_dist(0.0, 1.0);
        //std::vector<int> change = {11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 108, 109, 110, 111, 112, 113};
        std::vector<int> change = {0, 1, 4, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63, 67, 71, 75, 79, 83, 87, 91, 95, 99, 103};

        std::vector<json> genome;
        flatten(params, genome);

        for (int i = 0; i < change.size(); i++)
        {
            if (get_rand(0.0, 1.0) < mutation_rate)
            {
                genome[change[i]] = genome[change[i]].get<double>() + std::abs(dist(gen));
            }
        }

        genome[27] = genome[23];
        genome[35] = genome[31];
        genome[47] = genome[43];
        genome[55] = genome[51];
        genome[63] = genome[59];
        genome[71] = genome[67];
        genome[79] = genome[75];
        genome[87] = genome[83];
        genome[95] = genome[91];
        genome[103] = genome[99];

        json mutated_params;
        int idx = 0;
        unflatten(params, genome, idx, mutated_params);

        return mutated_params;

    }

    json write_base_params()
    {
        std::ifstream params_file("input/params.json");

        json params = json::parse(
                std::string(
                        (std::istreambuf_iterator<char>(params_file)),
                        (std::istreambuf_iterator<char>())
                ),
                nullptr,
                true,
                true
        );

        return params;
    }

    json write_random_params()
    {
        //std::vector<int> change = {11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 108, 109, 110, 111, 112, 113};
        std::vector<int> change = {0, 1, 4, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63, 67, 71, 75, 79, 83, 87, 91, 95, 99, 103};

        json base_params = write_base_params();
        std::vector<json> base_genome;
        flatten(base_params, base_genome);

        double alpha = 531242.7639479245;
        base_genome[0] = get_rand(0.5 * alpha, 1.5 * alpha);

        double beta = 1.8663559499166489;
        base_genome[1] = get_rand(0.5 * beta, 1.5 * beta);

        double gamma = 1.1947655088094953;
        base_genome[4] = get_rand(0.5 * gamma, 1.5 * gamma);

        double AWA_AIY = -17.613256896037917;
        base_genome[11] = get_rand(1.5 * AWA_AIY, 0.5 * AWA_AIY);

        double AIY_AIY = 9.986898926639817;
        base_genome[15] = get_rand(0.5 * AIY_AIY, 1.5 * AIY_AIY);

        double AIY_RIA = -7.0634496816364285;
        base_genome[19] = get_rand(1.5 * AIY_RIA, 0.5 * AIY_RIA);

        double RIA_RMD = -15.734477641464135;
        base_genome[23] = get_rand(1.5 * RIA_RMD, 0.5 * RIA_RMD);
        base_genome[27] = base_genome[23];

        double SMD_RIA = 0.7060810365603275;
        base_genome[31] = get_rand(0.5 * SMD_RIA, 1.5 * SMD_RIA);
        base_genome[35] = base_genome[31];

        double RIA_RIA = -0.35546625206879734;
        base_genome[39] = get_rand(1.5 * RIA_RIA, 0.5 * RIA_RIA);

        double RIA_SMD = 15.15488076456896;
        base_genome[43] = get_rand(0.5 * RIA_SMD, 1.5 * RIA_SMD);
        base_genome[47] = base_genome[43];

        double RMD_RIA = 5.927444644786412;
        base_genome[51] = get_rand(0.5 * RMD_RIA, 1.5 * RMD_RIA);
        base_genome[55] = base_genome[51];

        double SMD_SMD = -14.9121;
        base_genome[59] = get_rand(1.5 * SMD_SMD, 0.5 * SMD_SMD);
        base_genome[63] = base_genome[59];

        double RMD_RMD = 6.62512;
        base_genome[67] = get_rand(0.5 * RMD_RMD, 1.5 * RMD_RMD);
        base_genome[71] = base_genome[67];

        double SMD_SMDx = -11.2755;
        base_genome[75] = get_rand(1.5 * SMD_SMDx, 0.5 * SMD_SMDx);
        base_genome[79] = base_genome[75];

        double SMD_RMD = 14.9933;
        base_genome[83] = get_rand(0.5 * SMD_RMD, 1.5 * SMD_RMD);
        base_genome[87] = base_genome[83];

        double RMD_RMDx = -11.6075;
        base_genome[91] = get_rand(1.5 * RMD_RMDx, 0.5 * RMD_RMDx);
        base_genome[95] = base_genome[91];

        double SMDx_RMDx = 0.0199558;
        base_genome[99] = get_rand(0.5 * SMDx_RMDx, 1.5 * SMDx_RMDx);
        base_genome[103] = base_genome[99];

        //base_genome[108] = get_rand(0.0034, 1);
        //base_genome[109] = get_rand(-4, 5);
        //base_genome[110] = get_rand(0.004, 1);
        //base_genome[111] = get_rand(-10, -1);
        //base_genome[112] = get_rand(0.005, 1);
        //base_genome[113] = get_rand(-3, -8);

        int idx = 0;
        json random_params;
        unflatten(base_params, base_genome, idx, random_params);

        return random_params;

    }

};


#endif //GENETIC_ALGORITHM_H
