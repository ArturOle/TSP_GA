using Plots
using Random
using Logging
using Distributions
using BenchmarkTools


# The data structure to hold all relevant information about single idividual
mutable struct Indiv <:Number
    chromosome::Vector{Int64}
    fit::Float64
end

# Simple overlay over the vector to improve the code logics
mutable struct Generation <: Number
    individuals::Vector{Indiv}
end

Indiv(path) = Indiv(path, Inf)

function sub_distances(a::Indiv, b::Indiv)
    return a.dist - b.dist
end

function sub_distances(a::Number, b::Indiv)
    return a - b.dist
end

function output_function(cities, distances)
    distance = 0

    for i in 1:length(cities)-1
        distance += distances[cities[i], cities[i+1]]
    end
    distance += distances[cities[end], cities[1]]
    return distance
end

function evaluate_generation!(generation::Generation, distances)
    # For every individual
	for individual in generation.individuals
        # add evaluation at x(i) with a b and c of individual
        individual.fit = output_function(individual.chromosome, distances)
	end
    return generation
end

function evaluate_generation!(generation::Array, distances)
    # For every individual
	for individual in generation
        # add evaluation at x(i) with a b and c of individual
        individual.fit = output_function(individual.chromosome, distances)
	end
    return generation
end

# Initializes the first generation of individual according to the given population quantity
function initialize_population!(population, init_data, n=20, initialization=randperm)

    # We need to initialize the Vector of Indiv in form of the Generation
    append!(population, Generation(Vector{Indiv}(undef, n)))

    for i in 1:n
        population[1].individuals[i] = Indiv(
            initialization(init_data),              
            NaN
        )
    end
end

function sort_generation!(generation::Generation, desired_quantity::Int=200)
    sort!(generation.individuals, by=v->v.fit)[1:desired_quantity]
end

function sort_generation(generation::Generation)
    generation.individuals = sort(generation.individuals, by=v->v.fit)
    return generation
end

function sort_generation(generation::Generation, desired_quantity)
    generation.individuals = sort(generation.individuals, by=v->v.fit)[1:desired_quantity]
end

function rulette_selection(generation::Generation, desired_quantity::Int=200)
    result = Vector{Indiv}()
    gen_len = length(generation.individuals)
    res_len = 0
    sort_generation!(generation, desired_quantity)
    generation.individuals = generation.individuals[1:desired_quantity]
    best = generation.individuals[1].fit

    while true
        for individual in generation.individuals
            probability = best / individual.fit
            chance = rand(Uniform(0.0, 1.0))
            if probability > chance
                append!(result, individual)
                res_len += 1
            end
            if res_len >= gen_len
                return result[1:gen_len]
            end
        end
    end
end

function new_generation_evo(distances, population, selected, mutation_probability)
    # display(population[1])
    # Make crossover based on selected data and generate 5*population quantity of children
    offspring = crossover(length(population[1].individuals), selected)
    # Mutate all of children
    offspring = mutation(offspring, mutation_probability)
    # Mutate parents
    selected = mutation(selected, mutation_probability)
    # Evaluate new generation
    offspring = evaluate_generation!(offspring, distances)
    # Select the best from the population λ, γ and return
    return sort_generation(Generation(Vector{Indiv}(vcat(offspring, selected))))
end

function crossover(len_gen, selected)
    len_s = length(selected)
    offspring = []

    # We are generating 5*difference between the sizes of the base population 
    for i in 1:len_gen*5
        # Choosing the first parent randomly
        parent1 = rand(1:len_s)
        # Choosing the second parent randomly from population without parent1
        leftover = [r for r in 1:len_s if r!=parent1]
        if isempty(leftover)
            parent2=parent1
        else
            parent2 = rand(leftover)
        end
        # Creating Child
        child_1 = cross_two(selected[parent1], selected[parent2])
        append!(offspring, child_1)
        child_2 = cross_two(selected[parent2], selected[parent1])
        append!(offspring, child_2)
    end
    # Return the offspring
    return offspring
end

function cross_two(parent_first, parent_second)
    # Creating the child based on the prroperties of the parents
    
    nr_of_gens = length(parent_first.chromosome)

    child_chromosome = zeros(Int64, nr_of_gens)
    child_chromosome[1] = parent_first.chromosome[1]
    values = [child_chromosome[1]]
    positions = [1]
    val = parent_second.chromosome[2]
    pos = findfirst(x -> x==val, parent_second.chromosome)

    leftover_cities = []
    empty_positions = []

    if val in values || pos in positions
        empty_positions = [x[1] for x in enumerate(child_chromosome) if x[2]==false]
        if isempty(empty_positions)
            return Indiv(child_chromosome, NaN)
        end
        leftover_cities = [y for y in parent_first.chromosome if y ∉ values]
        child_chromosome[empty_positions] .= leftover_cities
    else
        append!(values, val)
        append!(positions, pos)
        child_chromosome[positions[end]] = val
    end

    for i in 2:nr_of_gens

        val = parent_first.chromosome[i]
        pos = findfirst(x -> x==val, parent_first.chromosome)

        if val in values || pos in positions
            empty_positions = [x[1] for x in enumerate(child_chromosome) if x[2]==false]
            if isempty(empty_positions)
                break
            end
            leftover_cities = [y for y in parent_second.chromosome if y ∉ values]
            child_chromosome[empty_positions] .= leftover_cities
        else
            append!(values, val)
            append!(positions, pos)
            child_chromosome[positions[end]] = val
        end

        val = parent_second.chromosome[i]
        pos = findfirst(x -> x==val, parent_second.chromosome)

        if val in values || pos in positions
            empty_positions = [x[1] for x in enumerate(child_chromosome) if x[2]==false]
            if isempty(empty_positions)
                break
            end
            leftover_cities = [y for y in parent_second.chromosome if y ∉ values]
            child_chromosome[empty_positions] .= leftover_cities
        else
            append!(values, val)
            append!(positions, findfirst(x -> x==val, parent_second.chromosome))
            child_chromosome[positions[end]] = val
        end
    end

    return Indiv(child_chromosome, NaN)
end


function mutation(offspring, mutation_probability)

    nr_of_genes = length(offspring[1].chromosome)
    nr_of_offsprings = length(offspring)

    for i in 1:nr_of_offsprings
        if mutation_probability > rand(Uniform(0.0, 1.0))
            pair = rand(1:nr_of_genes, 2)
            swap!(offspring[i].chromosome[pair[1]], offspring[i].chromosome[pair[2]])
        end
    end

    return offspring
end

function swap!(a, b)
    a, b = b, a
end

function EvolutionAlgorithm(
    data,
    population_quantity::Int=100,
    epsilon=0.00001,
    mutation_probability=0.2
    )

    N = length(data[1])
    generation = 1
    population = Vector{Generation}()
    
    distances = zeros(Float64, (N, N))

    for i in 1:N
        for j in i:N
            dist = sqrt((data[1][i] - data[1][j])^2 + (data[2][i] - data[2][j])^2)
            distances[i,j] = dist
            distances[j,i] = dist
        end
    end

    initialize_population!(population, N, population_quantity)
    evaluate_generation!(population[1], distances)

    while generation < population_quantity

        selected = rulette_selection(population[generation], population_quantity)
        next_generation = new_generation_evo(distances, population[generation], selected, mutation_probability)
        generation += 1
        append!(population, next_generation)
        new_best = population[generation].individuals[1].fit
        best = mean([x.fit for x in population[generation-1].individuals[1:floor(Int, population_quantity/10)]])

        if abs(new_best - best) > epsilon
            best = new_best
        else
            evaluate_generation!(
                population[end],
                distances
            )
            break
        end

        evaluate_generation!(
            population[end],
            distances
        )
        
    end


    return [select_parents(population, generation), generation]
end

function select_parents(population, generation=1)
    population[generation].individuals = sort(population[generation].individuals, by=v -> v.fit)
    number = Int(floor(length(population[1].individuals)/10))
	return population[generation].individuals[1:number]
end

function visualize_graph(data, best)
	points_coords_x = data[1]
	points_coords_y = data[2]

	plot(legend=false)

	for index in 2:length(data[1])
		first_point = (data[1][index], data[2][index])
		for jndex in 1:index-1
			second_point = (data[1][jndex], data[2][jndex])
			x = [first_point[1], second_point[1]]
			y = [first_point[2], second_point[2]]
			path_weight = 1

			plot!(x, y, lw=path_weight*1, color="lightblue")
		end
	end

    for i in 1:length(data[1])-1  #(i, path) in enumerate(zip(best))
        plot!(
            [data[1][best.chromosome[i]], data[1][best.chromosome[i+1]]],
            [data[2][best.chromosome[i]], data[2][best.chromosome[i+1]]],
            lw=3,
            color="#73C6B6"
        )   
    end
    plot!(
        [data[1][best.chromosome[end]], data[1][best.chromosome[1]]],
        [data[2][best.chromosome[end]], data[2][best.chromosome[1]]],
        lw=3,
        color="#73C6B6"
    ) 
	scatter!(points_coords_x, points_coords_y, color="#73C6B6", series_annotations=text.(1:length(data[1]), :bottom))
end


# x = [3 2 12 7  9  3 16 11 9 2]
# y = [1 4 2 4.5 9 1.5 11 8 10 7]
x = [0, 3, 6, 7, 15, 12, 14, 9, 7, 0]
y = [1, 4 ,5, 3, 0, 4, 10, 6, 9, 10]
data = [x,y]
@benchmark EvolutionAlgorithm(data)
# display(best[1][1])
# visualize_graph(data, best[1][1])
