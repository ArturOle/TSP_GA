using Test
include("tsp_ga.jl")


@testset "test cross_two" begin
    parent_a = Indiv([1,2,3,4,5], NaN)
    parent_b = Indiv([5,4,3,2,1], NaN)

    child = cross_two(parent_a, parent_b)
    @test child.chromosome == [1, 4, 5, 3, 2]

    parent_a = Indiv([1,2,3,4,5], NaN)
    parent_b = Indiv([1,2,3,4,5], NaN)

    child = cross_two(parent_a, parent_b)
    @test child.chromosome == [1,2,3,4,5]

    parent_a = Indiv([1,4,3,2,5], NaN)
    parent_b = Indiv([5,1,3,4,2], NaN)

    child = cross_two(parent_a, parent_b)
    @test child.chromosome == [1, 4, 3, 2, 5]

    parent_a = Indiv([1,5,3,2,4], NaN)
    parent_b = Indiv([5,1,3,4,2], NaN)

    child = cross_two(parent_a, parent_b)
    @test child.chromosome == [1, 5, 3, 2, 4]
end


@testset "test cross_two_cycle" begin
    parent_a = Indiv([1,2,3,4,5], NaN)
    parent_b = Indiv([5,4,3,2,1], NaN)

    child = cross_two_cycle(parent_a, parent_b)
    @test child.chromosome == [1, 4, 3, 2, 5]

    parent_a = Indiv([1,2,3,4,5], NaN)
    parent_b = Indiv([1,2,3,4,5], NaN)

    child = cross_two_cycle(parent_a, parent_b)
    @test child.chromosome == [1, 2, 3, 4, 5]

    parent_a = Indiv([1,4,3,2,5], NaN)
    parent_b = Indiv([5,1,3,4,2], NaN)

    child = cross_two_cycle(parent_a, parent_b)
    @test child.chromosome == [1, 3, 4, 2, 5]

    parent_a = Indiv([1,5,3,2,4], NaN)
    parent_b = Indiv([5,1,3,4,2], NaN)

    child = cross_two_cycle(parent_a, parent_b)
    @test child.chromosome == [1, 5, 3, 4, 2]
end
