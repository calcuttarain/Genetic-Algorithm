# Genetic Algorithm for Function Optimization

This repository contains Python code implementing a genetic algorithm for function optimization. The algorithm is an academic example of the genetic algorithms usage designed to find the maximum value of a quadratic function within a given domain.

## Features

- **Genetic Algorithm**: Utilizes a genetic algorithm approach to evolve a population of chromosomes towards the optimal solution.
- **Function Optimization**: Seeks to find the maximum value of a quadratic function within a specified domain.
- **Selection, Crossover, Mutation**: Implements selection, crossover, and mutation operations to iteratively improve the population.
- **Fitness Evaluation**: Evaluates the fitness of each chromosome based on the quadratic function.
- **Visualization**: Provides visualization of the evolution process, showing the mean fitness and the maximum fitness achieved over generations.

## Requirements

- Python 3.x
- Matplotlib

## Usage

1. Clone the repository:

    ```
    git clone https://github.com/yourusername/genetic-algorithm.git
    ```

2. Navigate to the project directory:

    ```
    cd genetic-algorithm
    ```

3. Ensure your input file (`input.txt`) containing algorithm parameters and domain specifications is in the root directory. The input file should have the following format:
    - Number of chromosomes
    - Domain of the function (closed interval endpoints)
    - Coefficients of the quadratic function
    - Precision for discretizing the interval
    - Crossover probability
    - Mutation probability
    - Number of generations
    
    Refer to the provided `input.txt` file for an example.

4. Run the Python script:

    ```
    python genetic_algorithm.py
    ```

5. Check the output file (`output.txt`) for the results, including initial and final populations, as well as statistics on the evolution process.

## Output

The algorithm generates an output file (`output.txt`) containing:

- Initial population
- Selection probabilities for the initial population
- Chromosomes selected for crossover
- Chromosomes after crossover
- Chromosomes selected for mutation
- Chromosomes after mutation
- Statistics on mean fitness and maximum fitness over generations

