#include "genetic_algorithm.h"

int main()
{
    GeneticAlgorithm genetic_algorithm;

    std::setprecision(10);

    double fit;
    double velocity_weight = 0.0; double chemosensation_weight = 1.0;
    std::array<json, individuals> gen_params;
    std::vector<double> fitness_avg(individuals, 0);
    std::vector<double> proximity_fitness(individuals, 0);
    std::vector<double> velocity_fitness(individuals, 0);
    std::vector<double> chemosensation_fitness(individuals, 0);
    std::array<json, individuals / 2> final_params;
    json params;
    //std::vector<double> foodpos_x = {-0.01, -0.01 * (std::sqrt(2) / 2), 0.0, 0.01 * (std::sqrt(2) / 2), 0.01};
    //std::vector<double> foodpos_y = {0.0, 0.01 * (std::sqrt(2) / 2), 0.01, 0.01 * (std::sqrt(2) / 2), 0.0};
    std::vector<double> foodpos_x = {-0.0042426, -0.0042426 * (std::sqrt(2) / 2), 0.0, 0.0042426 * (std::sqrt(2) / 2), 0.0042426};
    std::vector<double> foodpos_y = {0.0, 0.0042426 * (std::sqrt(2) / 2), 0.0042426, 0.0042426 * (std::sqrt(2) / 2), 0.0};

    std::string output_file = "log/0.txt";
    std::ofstream output(output_file, std::ios::app);

    std::string fit_file = "log/fitness.txt";
    std::ofstream fit_f(fit_file, std::ios::app);

    std::string diversity_file = "log/diversity.txt";
    std::ofstream diversity_f(diversity_file, std::ios::app);

    // Run first generation.

    for (int i = 0; i < individuals; i++)
    {

        // Hybrid initialization.
        if (i == 0)
        {
            params = genetic_algorithm.write_random_params();
        } else
        {
            // Randomize parameters for individuals.
            params = genetic_algorithm.write_random_params();
        }

        // Run two simulations per individual at two different initial angles.
        for (int j = 0; j < 5; j++)
        {
            params["ChemoReceptors"]["foodPos"]["x"] = foodpos_x[j];
            params["ChemoReceptors"]["foodPos"]["y"] = foodpos_y[j];
            params["ChemoReceptors"]["kappa"] = genetic_algorithm.get_rand(-1.0, -0.1);

            // Set random seed.
            json simulation_params = params["simulation"];
            long seed = genetic_algorithm.set_seed(simulation_params);
            RandomState rs;
            rs.SetRandomSeed(seed);

            // Load collision objects.
            std::vector<CollisionObject> collObjs = load_objects(params["simulation"]["coll"]);

            // Initialize worm object.
            InitializeBodyConstants();
            Worm w(params);

            // Run simulation.
            double distancetravelled;
            double fitness_concentration;
            std::tie(distancetravelled,
                     fitness_concentration) = genetic_algorithm.EvaluationFunction(w,
                                                                                   rs,
                                                                                   params["simulation"]["angle"],
                                                                                   collObjs,
                                                                                   params["simulation"]["t_food_start"].get<double>(),
                                                                                   output_file,
                                                                                   foodpos_x[j],
                                                                                   foodpos_y[j]);

            // Get fitness.
            double velocity, chemosensation;
            std::tie(velocity,
                     chemosensation) = genetic_algorithm.fitness(distancetravelled,
                                                                 fitness_concentration);

            // Combine fitnesses for both simulations.
            fitness_avg[i] += 0.2 * (velocity_weight * velocity + chemosensation_weight * chemosensation);
            std::cout << "Worm " << i << "/" << individuals << "." << j << ": " << "Fitness = " << fitness_avg[i]
                      << " = " << velocity_weight << " * " << velocity
                      << " + " << chemosensation_weight << " * " << chemosensation << std::endl;
            gen_params[i] = params;
        }
    }

    // Select fittest individuals.
    std::array<json, individuals / 2> selected_params;
    for (int i = 0; i < individuals / 2; i++)
    {
        json parent;
        // Elitism to always select closest proximity.
        if (i < elitism_size)
        {
            std::vector<int> fit_vector = genetic_algorithm.get_fittest(chemosensation_fitness);
            std::cout << "Elitism " << i << "..." << std::endl;
            parent = genetic_algorithm.elitism_select(gen_params, chemosensation_fitness, i);
        } else
        {
            std::cout << "Fitness proportional selection " << i << "..." << std::endl;
            parent = genetic_algorithm.fitness_proportional_select(gen_params,
                                                                   fitness_avg);
        }
        selected_params[i] = parent;
    }

    final_params = selected_params;

    // Crossover.
    std::cout << "Crossover..." << std::endl;
    for (int i = 0; i < individuals; i++)
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> distrib(0, (individuals/2)-1);

        int rand_1 = distrib(gen);
        int rand_2 = distrib(gen);

        json params_1 = selected_params[rand_1];
        json params_2 = selected_params[rand_2];
        crossed[i] = genetic_algorithm.crossover(params_1,
                                                 params_2,
                                                 crossover_rate);
    }

    // Mutation.
    std::cout << "Mutation..." << std::endl;
    for (int i = 0; i < individuals; i++)
    {
        mutated[i] = genetic_algorithm.mutate(crossed[i],
                                              mutation_rate,
                                              0.01);
    }

    // Get diversity measure.
    double distance = genetic_algorithm.hamming_distance(gen_params);
    std::cout << "Hamming distance: " << distance << std::endl;
    output << "\nHamming distance: " << std::fixed << std::setprecision(10) << distance << std::endl;
    diversity_f << distance << std::endl;

    std::vector<double> fitness_to_print;
    double fitness_print = 0.0;
    std::vector<int> fitness_indices = genetic_algorithm.get_fittest(fitness_avg);
    for (int index: fitness_indices)
    {
        double fitness_value = fitness_avg[index];
        fitness_to_print.push_back(fitness_value);
    }
    for (int i = 0; i < individuals / 2; i++)
    {
        fitness_print += fitness_to_print[i];
    }
    fitness_print /= (individuals/2);

    std::cout << "Generation 0 fitness: " << fitness_print << std::endl;
    output << "\nFitness: " << fitness_print << std::endl;

    fit_f << fitness_print << std::endl;

    output << "\n========================" << std::endl;
    output << "--generation parameters" << std::endl;
    output << "========================\n" << std::endl;

    for (int i = 0; i < individuals; i++)
    {
        for (const auto& item : gen_params[i].items())
        {
            output << std::fixed << std::setprecision(10) << item.key() << ": " << item.value() << std::endl;
        }
    }

    output << "\n===================" << std::endl;
    output << "--sorted fitnesses" << std::endl;
    output << "===================\n" << std::endl;

    for (int value: genetic_algorithm.get_fittest(fitness_avg))
    {
        output << std::fixed << std::setprecision(10) << fitness_avg[value] << std::endl;
    }

    output << "\n=====================" << std::endl;
    output << "--selected parameters" << std::endl;
    output << "=====================\n" << std::endl;

    for (int i = 0; i < individuals/2; i++)
    {
        for (const auto& item : selected_params[i].items())
        {
            output << std::fixed << std::setprecision(10) << item.key() << " " << item.value() << std::endl;
        }
    }

    output << "\n====================" << std::endl;
    output << "--crossed parameters" << std::endl;
    output << "====================\n" << std::endl;

    for (int i = 0; i < individuals; i++)
    {
        for (const auto& item : crossed[i].items())
        {
            output << std::fixed << std::setprecision(10) << item.key() << " " << item.value() << std::endl;
        }
    }

    output << "\n====================" << std::endl;
    output << "--mutated parameters" << std::endl;
    output << "====================\n" << std::endl;

    for (int i = 0; i < individuals; i++)
    {
        for (const auto& item : mutated[i].items())
        {
            output << std::fixed << std::setprecision(10) << item.key() << " " << item.value() << std::endl;
        }
    }

    double fittest = std::accumulate(fitness_avg.begin(),
                                     fitness_avg.end(), 0.0)
                     / fitness_avg.size();


    // Continue algorithm for remaining generations.
    for (int gen = 1; gen < gen_count; gen++)
    {
        std::array<json, individuals> gen_params;
        std::array<json, individuals / 2> selected_params_next;
        std::vector<double> fitness_avg_next(individuals, 0);
        std::vector<double> velocity_fitness_next(individuals, 0);
        std::vector<double> proximity_fitness_next(individuals, 0);
        std::vector<double> chemosensation_fitness_next(individuals, 0);
        std::vector<int> random_index;

        std::string output_file = "log/" + std::to_string(gen) + ".txt";
        std::ofstream output(output_file,std::ios::app);

        std::cout << "Generation " << gen << "/" << gen_count << "..." << std::endl;

        // Run simulation for each individual.
        for (int i = 0; i < individuals; i++)
        {

            // Run two simulations per individual with two different starting angles.
            for (int j = 0; j < 5; j++)
            {

                if (std::find(random_index.begin(), random_index.end(), i) != random_index.end())
                {
                    std::cout << "Randomizing " << i << "..." << std::endl;
                    params = genetic_algorithm.write_random_params();
                } else
                {
                    // Get parameters from previous generation.
                    params = mutated[i];
                }

                params["ChemoReceptors"]["foodPos"]["x"] = foodpos_x[j];
                params["ChemoReceptors"]["foodPos"]["y"] = foodpos_y[j];
                params["ChemoReceptors"]["kappa"] = genetic_algorithm.get_rand(-1.0, -0.1);

                json simulation_params = params["simulation"];

                // Set random seed.
                long seed = genetic_algorithm.set_seed(simulation_params);
                RandomState rs;
                rs.SetRandomSeed(seed);

                // Initialize worm.
                std::vector<CollisionObject> collObjs = load_objects(params["simulation"]["coll"]);
                InitializeBodyConstants();
                Worm w(params);

                // Store generation parameters.
                gen_params[i] = params;

                // Run simulation.
                double distancetravelled;
                double fitness_concentration;
                std::tie(distancetravelled,
                         fitness_concentration) =  genetic_algorithm.EvaluationFunction(w,
                                                                                        rs,
                                                                                        params["simulation"]["angle"],
                                                                                        collObjs,
                                                                                        params["simulation"]["t_food_start"].get<double>(),
                                                                                        output_file,
                                                                                        foodpos_x[j],
                                                                                        foodpos_y[j]);

                // Get fitness.
                double velocity, proximity, chemosensation;
                std::tie(velocity,
                         chemosensation) = genetic_algorithm.fitness(distancetravelled,
                                                                     fitness_concentration);

                // Combine fitnesses for both simulations per individual.
                fitness_avg_next[i] += 0.2 * (velocity_weight * velocity + chemosensation_weight * chemosensation);
                std::cout << "\r" << "Worm " << i << "/" << individuals << "." << j << ": " << "Fitness = " << fitness_avg_next[i]
                          << " = " << velocity_weight << " * " << velocity
                          << " + " << chemosensation_weight << " * " << chemosensation << std::endl;
            }

        }

        double distance_next = genetic_algorithm.hamming_distance(gen_params);
        std::cout << "Hamming distance:" << distance_next << std::endl;
        diversity_f << distance_next << std::endl;

        double this_fitness = std::accumulate(fitness_avg_next.begin(),
                                              fitness_avg_next.end(), 0.0)
                              / fitness_avg_next.size();

        // Adaptive crossover, mutation, and selection based on fitness and convergence.

        double rounded_this_fitness = std::round(this_fitness * 10000.0) / 10000.0;
        double rounded_fittest = std::round(fittest * 10000.0) / 10000.0;

        //if (this_fitness < fittest && distance_next < 0.2)
        //{
        //    std::cout << "Adjusting..." << std::endl;
        //    crossover_rate = 0.9;
        //    mutation_rate = 0.01;
        //    random_index = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
        //    elitism_size = std::round(individuals / 20.0);
        //} else
        //{
        //    fittest = this_fitness;
        //    crossover_rate = 0.6;
        //    mutation_rate = 0.005;
        //    random_index = {0};
        //    elitism_size = individuals / 5;
        // }

        // Select fittest individuals.
        for (int i = 0; i < individuals / 2; i++)
        {
            json parent;
            // Elitism based on proximity.
            if (i < elitism_size)
            {
                std::vector<int> fit_next = genetic_algorithm.get_fittest(chemosensation_fitness_next);
                std::cout << "Elitism " << i << "..." << std::endl;
                parent = genetic_algorithm.elitism_select(gen_params, chemosensation_fitness_next, i);
            } else
            {
                std::cout << "Fitness proportional selection " << i << "..." << std::endl;
                parent = genetic_algorithm.fitness_proportional_select(gen_params,
                                                                       fitness_avg_next);
            }
            selected_params_next[i] = parent;
        }

        // Crossover.
        std::cout << "Crossover..." << std::endl;
        for (int i = 0; i < individuals; i++)
        {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<int> distrib(0, (individuals/2)-1);

            int rand_1 = distrib(gen);
            int rand_2 = distrib(gen);

            json params_1 = selected_params_next[rand_1];
            json params_2 = selected_params_next[rand_2];
            crossed[i] = genetic_algorithm.crossover(params_1,
                                                     params_2,
                                                     crossover_rate);
        }

        // Mutation.
        std::cout << "Mutation..." << std::endl;
        for (int i = 0; i < individuals; i++)
        {
            mutated[i] = genetic_algorithm.mutate(crossed[i],
                                                  mutation_rate,
                                                  0.01);
        }

        output << "\nHamming distance: " << std::fixed << std::setprecision(10) << distance_next << std::endl;

        std::vector<double> fitness_to_print_next;
        double fitness_print_next = 0.0;
        std::vector<int> fitness_indices_next = genetic_algorithm.get_fittest(fitness_avg_next);
        for (int index: fitness_indices_next)
        {
            double fitness_value = fitness_avg_next[index];
            fitness_to_print_next.push_back(fitness_value);
        }
        for (int i = 0; i < individuals / 2; i++)
        {
            fitness_print_next += fitness_to_print_next[i];
        }
        fitness_print_next /= (individuals/2);

        std::cout << "Generation " << gen << " Fitness: " << fitness_print_next << std::endl;
        output << "\nFitness: " << fitness_print_next << std::endl;
        fit_f << fitness_print_next << std::endl;

        output << "\n=======================" << std::endl;
        output << "--generation parameters" << std::endl;
        output << "=======================\n" << std::endl;

        for (int i = 0; i < individuals; i++)
        {
            for (const auto& item : gen_params[i].items())
            {
                output << std::fixed << std::setprecision(10) << item.key() << ": " << item.value() << std::endl;
            }
        }

        output << "\n==================" << std::endl;
        output << "--sorted fitnesses" << std::endl;
        output << "==================\n" << std::endl;

        for (int value: genetic_algorithm.get_fittest(fitness_avg_next))
        {
            output << std::fixed << std::setprecision(10) << fitness_avg_next[value] << std::endl;
        }

        output << "\n=====================" << std::endl;
        output << "--selected parameters" << std::endl;
        output << "=====================\n" << std::endl;

        for (int i = 0; i < individuals/2; i++)
        {
            for (const auto& item : selected_params_next[i].items())
            {
                output << std::fixed << std::setprecision(10) << item.key() << " " << item.value() << std::endl;
            }
        }

        output << "\n====================" << std::endl;
        output << "--crossed parameters" << std::endl;
        output << "====================\n" << std::endl;

        for (int i = 0; i < individuals; i++)
        {
            for (const auto& item : crossed[i].items())
            {
                output << std::fixed << std::setprecision(10) << item.key() << " " << item.value() << std::endl;
            }
        }

        output << "\n=====================" << std::endl;
        output << "--mutated parameters" << std::endl;
        output << "====================\n" << std::endl;

        for (int i = 0; i < individuals; i++)
        {
            for (const auto& item : mutated[i].items())
            {
                output << std::fixed << std::setprecision(10) << item.key() << " " << item.value() << std::endl;
            }
        }

        final_params = selected_params_next;

        std::string fittest_params = "log/params_" + std::to_string(gen) + ".json";
        std::ofstream fittest_output(fittest_params);
        fittest_output << selected_params_next[0].dump(4);
    }

    // Save parameters from final generation.
    for (int i = 0; i < individuals / 2; i++)
    {
        //std::cout << "Final Params " << i << ":\n";
        //for (const auto& item : final_params[i].items()) {
        //    std::cout << item.key() << ": " << item.value() << "\n";
        //}
        std::string final_params_file = "log/final_params_" + std::to_string(i) + ".json";
        std::ofstream final_output(final_params_file);
        final_output << final_params[i].dump(4);

    }

    return 0;
}

