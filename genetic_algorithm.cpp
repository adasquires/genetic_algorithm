#include "genetic_algorithm.h"
#include "params.cpp"

int main()
{

    double fit;
    std::vector<double> fitness_avg(individuals);
    std::vector<double> fitness_avg_next(individuals);

    for (int i = 0; i < 2; i ++)
    {
        for (int index = 0; index < individuals; index++)
        {

            std::string str = "log/log-0-" + std::to_string(index) + ".txt";
            std::ofstream output(str.c_str());
            freopen(str.c_str(), "w", stdout);

            write_random_params();
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

            json simulation_params = params["simulation"];
            long seed = set_seed(simulation_params);
            RandomState rs;
            rs.SetRandomSeed(seed);

            std::vector<CollisionObject> collObjs = load_objects(params["simulation"]["coll"]);
            InitializeBodyConstants();
            Worm w(params);

            gen_params[index] = parse_params(params);

            if (i == 0)
            {
                double angle = params["simulation"]["angle"].get<double>()
            }
            else
            {
                double angle = -(params["simulation"]["angle"].get<double>())
            }

            std::pair<double, double> headpos = EvaluationFunction(
                w,
                rs,
                angle,
                collObjs,
                params["simulation"]["output"].get<std::string>(),
                params["simulation"]["t_food_start"].get<double>(),
                std::string("0-" + std::to_string(index))
            );

            double fit = fitness(headpos.first, headpos.second);
            fitnesses[index].push_back(fit);
            fitness_avg[index] += fit;
            std::ofstream f("log/fitness.txt", std::ios::app);
            f << ": 0." << index << "\n";

            for (const pair<const std::string, double>& value : gen_params[index])
            {
                output << value.first << ": " << value.second << std::endl;
            }

            for (const double& value : fitnesses[index])
            {
                output << value << std::endl;
            }
        }
    }

    selected_params = select(gen_params, fitness_avg);

    std::ofstream selected_file("selected.txt");
    if (selected_file.is_open())
    {
        for (int i = 0; i < individuals/2; i++)
        {
            for (const auto& pair : selected_params[i])
            {
              selected_file << pair.first << " " << pair.second << std::endl;
            }
        }
    selected_file.close();
    }

    for (int index = 0; index < individuals/4; index++)
    {
        crossed[index] = crossover(selected_params[index],
                                   selected_params[index + individuals/4],
                                   param_names);
    }

    for (int index = individuals/4; index < individuals/2; index++)
    {
        crossed[index] = crossover(selected_params[index],
                                   selected_params[index - individuals/5],
                                   param_names);
    }

    std::ofstream crossed_file("crossover.txt");
    if (crossed_file.is_open())
    {
        for (int i = 0; i < individuals/2; i++)
        {
            for (const auto& pair : crossed[i])
            {
                crossed_file << pair.first << " " << pair.second << std::endl;
            }
        }
    crossed_file.close();
    }


    for (int index = 0; index < individuals; index++)
    {
        std::string output = "input/params-0-" + std::to_string(index) + "-.txt";
        mutated[index] = mutate(crossed[index], 0.05, output, 0.05, param_names);
    }

    std::ofstream mutated_file("mutate.txt");
    if (mutated_file.is_open())
    {
        for (int i = 0; i < individuals/2; i++)
        {
           for (const auto& pair : mutated[i])
           {
                mutated_file << pair.first << " " << pair.second << std::endl;
           }
        }
    mutated_file.close();
    }

    for (int gen = 1; gen < gen_count; gen++)
    {
        for (int i = 0; i < 2; i ++)
        {
            for (int index = 0; index < individuals; index++)
            {

                std::string str = "log/log-" + std::to_string(gen) + "-" + std::to_string(index) + ".txt";
                std::ofstream output(str.c_str());
                freopen(str.c_str(), "w", stdout);

                write_random_params(i);
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

                json simulation_params = params["simulation"];
                long seed = set_seed(simulation_params);
                RandomState rs;
                rs.SetRandomSeed(seed);

                std::vector<CollisionObject> collObjs = load_objects(params["simulation"]["coll"]);
                InitializeBodyConstants();
                Worm w(params);

                gen_params[index] = parse_params(params);

                if (i == 0)
                {
                    double angle = params["simulation"]["angle"].get<double>()
                }
                else
                {
                    double angle = -(params["simulation"]["angle"].get<double>())
                }

                std::pair<double, double> headpos = EvaluationFunction(
                    w,
                    rs,
                    angle,
                    collObjs,
                    params["simulation"]["output"].get<std::string>(),
                    params["simulation"]["t_food_start"].get<double>(),
                    std::string(std::to_string(gen) + "-" + std::to_string(index))
                );

                double fit = fitness(headpos.first, headpos.second);
                fitnesses[index].push_back(fit);
                fitness_avg_next[index] += fit;
                std::ofstream f("log/fitness.txt", std::ios::app);
                f << ": " << gen + 1 << "." << index << "\n";

                for (const pair<const std::string, double>& value : gen_params[index])
                {
                    output << value.first << ": " << value.second << std::endl;
                }

                for (const double& value : fitnesses[index])
                {
                    output << value << std::endl;
                }
            }
        }

        selected_params = select(gen_params, fitness_avg_next);

        std::ofstream selected_file("selected.txt");
        if (selected_file.is_open())
        {
            for (int i = 0; i < individuals/2; i++)
            {
                for (const auto& pair : selected_params[i])
                {
                    selected_file << pair.first << " " << pair.second << std::endl;
                }
            }
        selected_file.close();
        }

        for (int index = 0; index < individuals/2; index++)
        {
        crossed[index] = crossover(selected_params[index],
                                   selected_params[index + individuals/4],
                                   param_names);
        }

        for (int index = individuals/4; index < individuals/2; index++)
        {
             crossed[index] = crossover(selected_params[index],
                                        selected_params[index - individuals/5],
                                        param_names);
        }

        std::ofstream crossed_file("crossover.txt");
        if (crossed_file.is_open())
        {
            for (int i = 0; i < individuals/2; i++)
            {
                for (const auto& pair : crossed[i])
                {
                    crossed_file << pair.first << " " << pair.second << std::endl;
                }
            }
        crossed_file.close();
        }

        for (int index = 0; index < individuals; index++)
        {
            std::string output = "input/params-" + std::to_string(gen) + "-" + std::to_string(index) + "-" + ".txt";
            mutated[index] = mutate(crossed[index], 0.05, output, 0.05, param_names);
        }

        std::ofstream mutated_file("mutate.txt");
        if (mutated_file.is_open())
        {
            for (int i = 0; i < individuals/2; i++)
            {
                for (const auto& pair : mutated[i])
                {
                    mutated_file << pair.first << " " << pair.second << std::endl;
                }
            }
        mutated_file.close();
        }
    }

    for (int index = 0; index < individuals; index++)
    {
        std::ofstream output("final_params.txt", std::ios::app);
        output << gen_params[index] << std::endl;
        output << fitness_avg_next[index] << std::endl;
        output << "\n" << std::endl;
        output.close();
    }

    return 0;
};

