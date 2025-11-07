#include "genetic_algorithm.h"
#include "params.cpp"

int main()
{

    GeneticAlgorithm genetic_algorithm;

    std::setprecision(10);

    std::array<std::map<std::string, double>, individuals> gen_params;
    std::vector<double> fitness_avg(individuals, 0);
    std::vector<double> proximity_fitness(individuals, 0);
    std::vector<double> velocity_fitness(individuals, 0);
    std::vector<double> chemosensation_fitness(individuals, 0);

    std::string output_file = "log/0.txt";
    std::ofstream output(output_file, std::ios::app);

    // Run first generation. 
    for (int i = 0; index < individuals; index++)
    {
        // Hybrid initialization.
        if (i == 0)
        {
            write_base_params();
        } else
        {
            // Randomize parameters for individuals.
            write_random_params();
        }

        // Run two simulations per individual at two different initial angles. 
        for (int j = 0; j < 2; j++)
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

            if (j == 0)
            {
                angle = 30.0;
                params["simulation"]["angle"] = 30.0;
            } else
            {
                angle = 30.0;
                params["simulation"]["angle"] = -30.0;
            }

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

            // Store generation parameters.
            gen_params[i] = parse_params(params, "input/params.json");

            if (i == 0)
            {
                double angle = params["simulation"]["angle"].get<double>()
            }
            else
            {
                double angle = -(params["simulation"]["angle"].get<double>())
            }

            // Run simulation.
            double distancetravelled, head_x, head_y, close, closest_time;
            std::vector<double>fitness_concentration;
            std::tie(distancetravelled,
                     head_x,
                     head_y,
                     close,
                     closest_time,
                     fitness_concentration) = genetic_algorithm.EvaluationFunction(w,
                                                                                   rs,
                                                                                   angle,
                                                                                   collObjs,
                                                                                   params["simulation"]["t_food_start"].get<double>(),
                                                                                   output_file);

            // Get fitness.
            double velocity, proximity, chemosensation;
            std::tie(velocity,
                     proximity,
                     chemosensation) = genetic_algorithm.fitness(distancetravelled,
                                                                 head_x,
                                                                 head_y,
                                                                 close,
                                                                 closest_time,
                                                                 fitness_concentration,
                                                                 "input/params.json");

            double fit1 = 0.5 * velocity + 0.5 * chemosensation;

            chemosensation_fitness[i] += 0.5*chemosensation;
            proximity_fitness[i] += 0.5*proximity;
            velocity_fitness[i] += 0.5*velocity;

            // Combine fitnesses for both simulations.
            fitness_avg[i] += 0.5*velocity + 0.5*chemosensation;
        }
    }

    // Select fittest individuals.
    std::array<std::map<std::string, double>, individuals / 2> selected_params;

    // Tournament selection. 
    selected_params = genetic_algorithm.tournament_select(gen_params, velocity_fitness, tournament_size, elitism_size);

    // Crossover.
    for (int i = 0; i < individuals; i++)
    {
        std::random_device rd;
        std::mt19937 gen_rand(rd());
        std::uniform_int_distribution<int> distrib(0, (individuals/2)-1);
        std::uniform_real_distribution<double> real_dist(0.0, 1.0);
        std::uniform_int_distribution<int> cp_dist(0, param_names.size() - 1);

        int rand_1 = distrib(gen_rand);
        int rand_2 = distrib(gen_rand);

        const std::map<std::string, double>& params_1 = selected_params[rand_1];
        const std::map<std::string, double>& params_2 = selected_params[rand_2];
        crossed[i] = genetic_algorithm.crossover(params_1,
                                                params_2,
                                                 param_names,
                                                 crossover_rate);
    }

    // Mutation.
    std::array<std::map<std::string, double>, individuals> mutated;
    for (int i = 0; i < individuals; i++)
    {
        mutated[i] = genetic_algorithm.mutate(crossed[i],
                                              mutation_rate,
                                              0.01,
                                              param_names,
											  change_params,
                                              output_file);
    }

    // Get diversity measure.
    double distance = genetic_algorithm.hamming_distance(gen_params);
    std::cout << "Hamming distance:" << distance << std::endl;
    output << "\nHamming distance: " << std::fixed << std::setprecision(10) << distance << std::endl;

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
    fitness_print /= (individuals/2.0);

    std::cout << "Fitness: " << fitness_print << std::endl;
    output << "\nFitness: " << fitness_print << std::endl;

    output << "\n========================" << std::endl;
    output << "--generation parameters" << std::endl;
    output << "========================\n" << std::endl;

    for (int i = 0; i < individuals; i++)
    {
        for (const pair<const std::string, double>& value : gen_params[i])
        {
            output << std::fixed << std::setprecision(10) << value.first << ": " << value.second << std::endl;
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
        for (const auto& pair : selected_params[i])
        {
           output << std::fixed << std::setprecision(10) << pair.first << " " << pair.second << std::endl;
        }
    }

    output << "\n====================" << std::endl;
    output << "--crossed parameters" << std::endl;
    output << "====================\n" << std::endl;

    for (int i = 0; i < individuals; i++)
    {
       for (const auto& pair : crossed[i])
        {
            output << std::fixed << std::setprecision(10) << pair.first << " " << pair.second << std::endl;
        }
    }

    output << "\n====================" << std::endl;
    output << "--mutated parameters" << std::endl;
    output << "====================\n" << std::endl;

    for (int i = 0; i < individuals; i++)
    {
        for (const auto& pair : mutated[i])
        {
            output << std::fixed << std::setprecision(10) << pair.first << " " << pair.second << std::endl;
        }
    }

    double fittest = std::accumulate(fitness_avg.begin(),
                                     fitness_avg.end(), 0.0)
                                     / fitness_avg.size();

    std::array<std::map<std::string, double>, individuals / 2> final_params;
    std::array<std::map<std::string, double>, individuals> gen_params_next = gen_params;

    // Continue algorithm for remaining generations. 
    for (int generation = 1; generation < gen_count; generation++)
    {
        std::array<std::map<std::string, double>, individuals / 2> selected_params_next;
        std::vector<double> fitness_avg_next(individuals, 0);
        std::vector<double> velocity_fitness_next(individuals, 0);
        std::vector<double> proximity_fitness_next(individuals, 0);
        std::vector<double> chemosensation_fitness_next(individuals, 0);
        std::vector<int> random_index;

        std::string output_file_next = "log/" + std::to_string(generation) + ".txt";
        std::ofstream output_next(output_file_next,std::ios::app);

        // Run simulation for each individual. 
        for (int index = 0; index < individuals; index++)
        {
            // Run two simulations per individual with two different starting angles. 
            for (int j = 0; j < 2; j++)
            {
                if (std::find(random_index.begin(), random_index.end(), i) != random_index.end())
                {
                    std::cout << "Randomizing..." << std::endl;
                    write_random_params();
                } else
                {
                // Get parameters from previous generation.
                    write_params(mutated[i],
                                 "input/params.json");
                }
                
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

                 if (j == 0)
                {
                    angle = 30.0;
                    params["simulation"]["angle"] = 30.0;
                } else
                {
                    angle = 30.0;
                    params["simulation"]["angle"] = -30.0;
                }

                json simulation_params = params["simulation"];

                // Set random seed. 
                long seed = set_seed(simulation_params);
                RandomState rs;
                rs.SetRandomSeed(seed);

                // Initialize worm. 
                std::vector<CollisionObject> collObjs = load_objects(params["simulation"]["coll"]);
                InitializeBodyConstants();
                Worm w(params);

                // Store generation parameters. 
                gen_params_next[i] = parse_params(params, "input/params.json");

                // Run simulation.
                double distancetravelled, head_x, head_y, close, closest_time;
            	std::vector<double>fitness_concentration;
            	std::tie(distancetravelled,
                         head_x,
                         head_y,
                         close,
                         closest_time,
                         fitness_concentration) =  genetic_algorithm.EvaluationFunction(w,
                                                                                        rs,
                                                                                        angle,
                                                                                        collObjs,
                                                                                        params["simulation"]["t_food_start"].get<double>(),
                                                                                        output_file_next);

                // Get fitness.
                double velocity, proximity, chemosensation;
            	std::tie(velocity,
                         proximity,
                         chemosensation) = genetic_algorithm.fitness(distancetravelled,
                                                                     head_x,
                                                   	                 head_y,
                                                   	                 close,
                                                   	                 closest_time,
                                                                     fitness_concentration,
                                                                     "input/params.json");

                // Combine fitnesses for both simulations per individual.
                double fit = 0.5 * velocity + 0.4 * proximity + 0.1 * chemosensation;
                proximity_fitness_next[i] += 0.5*proximity;
                velocity_fitness_next[i] += velocity;
                chemosensation_fitness_next[i] += chemosensation;
            }
        }

        double average_velocity = 0;
        for (int i = 0; i < velocity_fitness_next.size() / 2; i++) {
            average_velocity += velocity_fitness_next[i];
        }
        average_velocity /= (velocity_fitness_next.size() / 2.0);

        double average_chemosensation = 0;
        for (int i = 0; i < chemosensation_fitness_next.size() / 2; i++) {
            average_chemosensation += chemosensation_fitness_next[i];
        }
        average_chemosensation /= (chemosensation_fitness_next.size() / 2.0);

        if (average_velocity >= 0.4) {
            std::cout << "Sufficient velocity fitness reached" << std::endl;
            for (int i = 0; i < individuals; i++) {
                fitness_avg_next[i] = 0.5*velocity_fitness_next[i] + 0.5*chemosensation_fitness_next[i];
            }
        } else {
            std::cout << "Insufficient velocity fitness" << std::endl;
            for (int i = 0; i < individuals; i++) {
                fitness_avg_next[i] = 0.5*velocity_fitness_next[i] + 0.5*chemosensation_fitness_next[i];
            }
        }

        double distance_next = genetic_algorithm.hamming_distance(gen_params_next);
        std::cout << "Hamming distance:" << distance_next << std::endl;

        double this_fitness = std::accumulate(fitness_avg_next.begin(),
                                              fitness_avg_next.end(), 0.0)
                                              / fitness_avg_next.size();

        // Adaptive crossover, mutation, and selection based on fitness and convergence.
        double rounded_this_fitness = std::round(this_fitness * 10000.0) / 10000.0;
        double rounded_fittest = std::round(fittest * 10000.0) / 10000.0;

        if (this_fitness < fittest && distance_next < 0.2)
        {
            std::cout << "Adjusting..." << std::endl;
            crossover_rate = 0.9;
            mutation_rate = 0.005;
            random_index = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
            random_index = {0};
            elitism_size = std::round(individuals / 10);
            tournament_size = individuals / 20.0;
        } else
        {
            fittest = this_fitness;
            crossover_rate = 0.7;
            mutation_rate = 0.001;
            random_index = {0};
            elitism_size = individuals / 3;
            tournament_size = individuals / 10.0;
        }

        // Select fittest individuals.
        std::vector<int> fit_vector_next = genetic_algorithm.get_fittest(velocity_fitness_next);

        selected_params_next = genetic_algorithm.tournament_select(gen_params_next, velocity_fitness_next, tournament_size, elitism_size);

        // Crossover.
        for (int i = 0; i < individuals; i++)
        {
            std::random_device rd;
            std::mt19937 gen_rand_next(rd());
            std::uniform_int_distribution<int> distrib(0, (individuals/2)-1);
            std::uniform_real_distribution<double> real_dist(0.0, 1.0);
            std::uniform_int_distribution<int> cp_dist(0, param_names.size() - 1);

            int rand_1 = distrib(gen_rand_next);
            int rand_2 = distrib(gen_rand_next);

            const std::map<std::string, double>& params_1 = selected_params_next[rand_1];
            const std::map<std::string, double>& params_2 = selected_params_next[rand_2];

        crossed[i] = genetic_algorithm.crossover(params_1,
                                                 params_2,
                                                 param_names,
                                                 crossover_rate);
        }

        // Mutation.
        std::array<std::map<std::string, double>, individuals> mutated;
        for (int i = 0; i < individuals; i++)
        {
            mutated[i] = genetic_algorithm.mutate(crossed[i],
                                                  mutation_rate,
                                                  0.01,
                                                  param_names,
                                                  change_params,
                                                  output_file_next);
        }

        output_next << "\nHamming distance: " << std::fixed << std::setprecision(10) << distance_next << std::endl;

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
        fitness_print_next /= (individuals/2.0);

        std::cout << "Fitness: " << fitness_print_next << std::endl;
        output_next << "\nFitness: " << fitness_print_next << std::endl;

        output_next << "\n=======================" << std::endl;
        output_next << "--generation parameters" << std::endl;
        output_next << "=======================\n" << std::endl;

        for (int i = 0; i < individuals; i++)
        {
            for (const pair<const std::string, double>& value : gen_params_next[i])
            {
                output_next << std::fixed << std::setprecision(10) << value.first << ": " << value.second << std::endl;
            }
        }

        output_next << "\n==================" << std::endl;
        output_next << "--sorted fitnesses" << std::endl;
        output_next << "==================\n" << std::endl;

        for (int value: genetic_algorithm.get_fittest(fitness_avg_next))
        {
            output_next << std::fixed << std::setprecision(10) << fitness_avg_next[value] << std::endl;
        }

        output_next << "\n=====================" << std::endl;
        output_next << "--selected parameters" << std::endl;
        output_next << "=====================\n" << std::endl;

        for (int i = 0; i < individuals/2; i++)
        {
            for (const auto& pair : selected_params_next[i])
            {
                output_next << std::fixed << std::setprecision(10) << pair.first << " " << pair.second << std::endl;
            }
        }

        output_next << "\n====================" << std::endl;
        output_next << "--crossed parameters" << std::endl;
        output_next << "====================\n" << std::endl;

        for (int i = 0; i < individuals; i++)
        {
            for (const auto& pair : crossed[i])
            {
               output_next << std::fixed << std::setprecision(10) << pair.first << " " << pair.second << std::endl;
            }
        }

        output_next << "\n=====================" << std::endl;
        output_next << "--mutated parameters" << std::endl;
        output_next << "====================\n" << std::endl;

        for (int i = 0; i < individuals; i++)
        {
            for (const auto& pair : mutated[i])
            {
                output_next << std::fixed << std::setprecision(10) << pair.first << " " << pair.second << std::endl;
            }
        }

        final_params = selected_params_next;

        int chemosensation_fittest = 0;
        for (int i = 0; i < individuals; i++) {
            if (chemosensation_fitness_next[i] >= chemosensation_fitness_next[chemosensation_fittest])
            {
                chemosensation_fittest = i;
            }
        }

        std::string fittest_params = "log/params_" + std::to_string(generation) + ".json";
        write_params(selected_params_next[0], fittest_params);
    }

    // Save parameters from final generation.
    for (int i = 0; i < individuals/2; i++)
    {
            std::cout << "Final Params " << i << ":\n";
            for (const auto& pair : final_params[i]) {
                std::cout << pair.first << ": " << pair.second << "\n";
            }
        std::string final_params_file = "log/final_params_" + std::to_string(i) + ".json";
        write_params(final_params[i],
                     final_params_file);
    }

    return 0;
}

