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

#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

using nlohmann::json;

WormBody b; RandomState rs;
double angle = 1.437;
double t_food_start = 1.0;

double random_double;

double alpha; double beta; double foodPos_x; double foodPos_y; double gamma; double kappa; double lambda;

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

double headpos; double foodpos; double distance;

double duration = 100.0;

const int gen_count = 300; const int individuals = 60;

std::array<std::vector<double>, individuals> fitnesses;
std::array<std::map<std::string, double>, individuals> gen_params;
std::array<std::map<std::string, double>, individuals/2> selected_params;
std::array<std::map<std::string, double>, individuals> crossed;
std::array<std::map<std::string, double>, individuals> mutated;

inline int set_seed(json & simulation_params)
{
    long seed = static_cast<long>(time(NULL));
    return seed;
}

std::pair<double, double> EvaluationFunction(Worm w,
                               RandomState &rs,
                               double angle,
                               std::vector<CollisionObject> & collObjs,
                               std::string output_dir,
                               double t_food_start,
                               std::string output_file)
{

  std::ofstream output(output_dir + "log-" + output_file + ".txt");
  output << "  > opening output files in " << output_dir.c_str() << std::endl;

  std::ofstream fitfile;
  fitfile.open(output_dir + "fitness.yaml");
  std::ofstream bodyfile;
  bodyfile.open(output_dir + "body.dat");
  std::ofstream actfile;
  actfile.open(output_dir + "act.dat");
  std::ofstream curvfile;
  curvfile.open(output_dir + "curv.dat");

  output << "  > initializing curvature measurement" << std::endl;
  TVector<double> curvature(1, N_curvs);
  TVector<double> antpostcurv(1, 2);
  antpostcurv.FillContents(0.0);

  FitnessCalc fcalc(w);

  output << "  > initializing worm state" << std::endl;
  w.InitializeState(rs, angle, collObjs);

  output << "  > starting time loop:" << std::endl;
  for (double t = 0.0; t <= DURATION; t += STEPSIZE)
  {
      if (t > t_food_start)
      {
          w.chemo_re.enabled = true;
      }

      if ( (t - (int) t < STEPSIZE))
      {
          output << "    >>  time: " << t << " / " << DURATION << std::endl;
          std::string pos_output = w.strprintf();
          std::replace(pos_output.begin(), pos_output.end(), '\n', '\t');
          output << pos_output << std::endl;
      }

      w.Step(STEPSIZE, 1);
      output << "step: " << t << std::endl;

      fcalc.update();
      output << "fcalc: " << t << std::endl;

      w.Curvature(curvature);
      output << "curvature: " << curvature << std::endl;
      curvfile << curvature << endl;

      w.DumpBodyState(bodyfile, skip);

      w.DumpActState(actfile, skip);

  }

  double head_x = w.CoMx();
  double head_y = w.CoMy();

  output << "\n\n  > finished time loop!\n" << std::endl;
  output << "  > closing files, saving to\n" << output_dir.c_str() << std::endl;

  bodyfile.close();
  actfile.close();
  curvfile.close();
  fitfile << fcalc.strprintf();

  output << "\n\n    > fitnesses:\n" << std::endl;
  std::string fcalc_output = fcalc.strprintf();
  std::replace(fcalc_output.begin(), fcalc_output.end(), '\n', '\t');
  output << fcalc_output.c_str() << std::endl;

  return std::make_pair(head_x, head_y);

}

double get_rand(double min, double max)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);
    double rand_num;
    do
    {
        rand_num = dis(gen);
    } while (rand_num == 0.0);
    return rand_num;
}

double fitness(double head_x, double head_y)
{

    std::ifstream f("input/params.json");

    if (!f)
    {
        std::cerr << "params.json not found" << std::endl;
    }

    json data;
    f >> data;

    foodPos_x = data["ChemoReceptors"]["foodPos"]["x"].get<double>();
    foodPos_y = data["ChemoReceptors"]["foodPos"]["y"].get<double>();

    std::ifstream h("log/pos.txt");

    std::string line;
    while (std::getline(h, line))
    {
        std::cout << line << std::endl;
    }

    double dx = head_x - foodPos_x;
    double dy = head_y - foodPos_y;

    double distance = std::sqrt(dx * dx + dy * dy);

    std::ofstream g("log/fitness.txt", std::ios::app);
    g << distance;

    return distance;

    };

std::vector<int> get_fittest(const std::vector<double>& arr)
{
    if (arr.size() <= individuals/2)
    {
        std::vector<int> indices(arr.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(),
                  indices.end(),
                      [&](int i, int j){ return arr[i] < arr[j]; });
        return indices;
    }

    std::vector<int> indices(arr.size());
    std::iota(indices.begin(), indices.end(), 0);

    std::partial_sort(indices.begin(),
                    indices.begin() + individuals/2,
                      indices.end(),
                    [&](int i, int j)
                      { return arr[i] < arr[j]; });

    indices.resize(individuals/2);

    std::ofstream get_fittest_file("fittest.txt");

    if (get_fittest_file.is_open())
        {
        for (const auto& pair : indices)
        {
            get_fittest_file << pair << std::endl;
        }
    get_fittest_file.close();
    }

    std::ofstream i("indices.txt");
    for (const auto& pair : indices)
    {
        i << pair << std::endl;
    }
    return indices;
    };

std::array<std::map<std::string, double>, individuals/2> select(const std::array<std::map<std::string, double>, individuals>& params,
                                                                const std::vector<double> fit)
{
    std::vector<int> f = get_fittest(fit);
    std::array<std::map<std::string, double>, individuals/2> fittest;
    for (int i = 0; i < individuals/2; i++)
    {
        fittest[i] = params[f[i]];
    }
    std::ofstream select_file("select.txt");
    if (select_file.is_open())
    {
        for (int i = 0; i < individuals/2; i++)
        {
            for (const auto& pair : fittest[i])
            {
                select_file << pair.first << " " << pair.second << std::endl;
            }
        }
        select_file.close();
    }
    return fittest;
    };

std::map<std::string, double> crossover(const std::map<std::string, double> params_1,
                                        const std::map<std::string, double> params_2,
                                        std::vector<std::string> param_names)
{
    std::map<std::string, double> child;
    int crossover_point = std::rand() % param_names.size();
    for (int i = 0; i < crossover_point; ++i)
    {
        std::string p = param_names[i];
        child[p] = params_1.find(p) -> second;
    }
    for (int i = crossover_point; i < params_1.size(); ++i) {
        std::string p = param_names[i];
        child[p] = params_2.find(p) -> second;
    }
    std::ofstream crossover_file("crossover.txt");

    if (crossover_file.is_open())
    {
        for (const auto& pair : child) {
            crossover_file << pair.first << " " << pair.second << std::endl;
        }
        crossover_file.close();
    }
    return child;
    };

std::map<std::string, double> mutate(std::map<std::string, double>& params,
                                    double mutation_rate, std::string output,
                                    double mutation_sigma, std::vector<std::string> param_names)
{
    std::map<std::string, double> mutated_geno;
    mutated_geno = params;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dist(0.0, mutation_sigma);
    std::uniform_real_distribution<> prob_dist(0.0, 1.0);
    for (int i = 0; i < param_names.size(); ++i)
    {
        if (get_rand(0.0, 1.0) < mutation_rate)
        {
            std::string p = param_names[i];
            mutated_geno[p] += dist(gen);
        }
    }

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
    mutated_geno["DB_DB"] = mutated_geno["VBA_VBA"] = mutated_geno["VBP_VBP"];
    mutated_geno["DD_DD"] = mutated_geno["VDA_VDA"] = mutated_geno["VDP_VDP"];
    mutated_geno["DB_DD"] = mutated_geno["VBA_VDA"] = mutated_geno["VBP_VDP"];
    mutated_geno["DB_VDP"] = mutated_geno["DB_VDP"];
    mutated_geno["VBA_DD"] = mutated_geno["VBP_DD"];
    mutated_geno["DD_VDA"] = mutated_geno["DD_VDP"];
    mutated_geno["VDA_VDP"] = mutated_geno["VDA_VDP"];
    mutated_geno["fwd_DB_DB"] = mutated_geno["fwd_VBP_VBA"];
    mutated_geno["fwd_DD_DD"] = mutated_geno["fwd_VDP_VDA"];
    mutated_geno["SMDD_RMDD_ele"] = std::abs(mutated_geno["SMDD_RMDD_ele"]);
    mutated_geno["SMDV_RMDV_ele"] = std::abs(mutated_geno["SMDV_RMDV_ele"]);
    mutated_geno["RMDV_RMDD_ele"] = std::abs(mutated_geno["RMDV_RMDD_ele"]);
    mutated_geno["DD_VDA_ele"] = std::abs(mutated_geno["DD_VDA_ele"]);
    mutated_geno["DD_VDP_ele"] = std::abs(mutated_geno["DD_VDP_ele"]);
    mutated_geno["VDA_VDP_ele"] = std::abs(mutated_geno["VDA_VDP_ele"]);
    mutated_geno["VBA_VBP_ele"] = std::abs(mutated_geno["VBA_VB_eleP"]);
    mutated_geno["fwd_DB_DB"] = std::abs(mutated_geno["fwd_DB_DB"]);
    mutated_geno["fwd_VBP_VBA"] = std::abs(mutated_geno["fwd_VBP_VBA"]);
    mutated_geno["fwd_DD_DD"] = std::abs(mutated_geno["fwd_DD_DD"]);
    mutated_geno["fwd_VDP_VDA"] = std::abs(mutated_geno["fwd_VDP_VDA"]);
    mutated_geno["fwd_VBP_DB"] = std::abs(mutated_geno["fwd_VBP_DB"]);

    std::ofstream mutate_file("mutate.txt");

    if (mutate_file.is_open())
    {
        for (const auto& pair : mutated_geno)
        {
            mutate_file << pair.first << " " << pair.second << std::endl;
        }
        mutate_file.close();
    }
    return mutated_geno;
  };

  #endif //GENETIC_ALGORITHM_H
