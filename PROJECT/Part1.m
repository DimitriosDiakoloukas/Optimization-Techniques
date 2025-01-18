clc;
clear;

function [weight] = GenerateInitialPopulation(V)
    weight = zeros(1, 17);
    function [val1, val2] = validSplit(total, limit1, limit2)
        while true
            val1 = rand * total;
            val2 = total - val1;
            if val1 < limit1 && val2 < limit2 && all([val1, val2] >= 0)
                break;
            end
        end
    end

    function values = validMultiSplit(total, limits)
        while true
            values = rand(1, numel(limits) - 1) .* total; 
            values(end + 1) = total - sum(values);       
            if all(values < limits) && all(values > 0)
                break;
            end
        end
    end

    while true
        weight(1:4) = [rand * 54.13, rand * 21.56, rand * 34.08, rand * 49.19];
        total = sum(weight(1:4));
        weight(1:4) = (weight(1:4) / total) * V;

        if all(weight(1:4) < [54.13, 21.56, 34.08, 49.19]) && all(weight(1:4) >= 0)
            break;
        end
    end

    [weight(5), weight(6)] = validSplit(weight(1), 33.03, 21.84);
    [weight(7), weight(8)] = validSplit(weight(2), 29.96, 24.87);
    [weight(9), weight(10)] = validSplit(weight(4), 47.24, 33.97);

    combinedFlow = weight(9) + weight(3) + weight(8);
    while true
        weight(11:13) = validMultiSplit(combinedFlow, [26.89, 32.76, 39.98]);
        weight(17) = weight(10) + weight(11);

        if weight(17) < 59.73
            break;
        end
    end

    combinedFlow2 = weight(6) + weight(7) + weight(13);
    while true
        weight(14) = rand * 37.12;
        weight(15) = combinedFlow2 - weight(14);
        weight(16) = weight(5) + weight(14);

        if all(weight([14, 15, 16]) < [37.12, 53.83, 61.65]) && all(weight([14, 15, 16]) > 0)
            break;
        end
    end
end


function [f] = FitnessFunc(a, c, currentPopulation, t, tolerance)
    currentPopulation = max(min(currentPopulation, c - tolerance), 0); 
    time = zeros(size(currentPopulation));
    total = zeros(1, size(currentPopulation, 1));
    time = t + a .* (currentPopulation ./ (1 - currentPopulation ./ c));
    total = sum(time, 2)'; 
    total = max(total, 1e-6);
    f = 1 ./ total';
end


function [cnt, I] = PerformRouletteSelection(fit)
    total = sum(fit);
    prob = fit / total;
    [prob_new, I] = sort(prob);
    q = cumsum(prob_new); 
    cnt = zeros(10, 1);
    randomValues = rand(10, 1);
    for k = 1:10
        idx = find(randomValues(k) < q, 1);
        cnt(idx) = cnt(idx) + 1;
    end
end


function [newPopulation] = GenerateNextPopulation(currentPopulation, fitnessProportions, selectedIndices)
    newPopulation = zeros(10, size(currentPopulation, 2));
    count = 1;
    for i = 1:length(fitnessProportions)
        if fitnessProportions(i) > 0
            numCopies = fitnessProportions(i); 
            for j = 1:numCopies
                newPopulation(count, :) = currentPopulation(selectedIndices(i), :);
                count = count + 1;
            end
        end
    end
end


function [weight_new] = ApplyMutation(weight, c, V, tolerance)
    weight_new = weight; 
    mu = 0;       
    sigma = 1;   
    numBranchesToChange = randi([2, 4]);
    branchesToMutate = randperm(4, numBranchesToChange);
    changes = (rand(1, numBranchesToChange) - 0.5) .* normpdf(weight(branchesToMutate), mu, sigma);
    for i = 1:numBranchesToChange
        weight_new(branchesToMutate(i)) = weight(branchesToMutate(i)) + changes(i);
    end

    totalChange = sum(changes); 
    remainingBranches = setdiff(1:4, branchesToMutate); 

    if ~isempty(remainingBranches)
        redistribution = rand(1, length(remainingBranches));
        redistribution = redistribution / sum(redistribution) * -totalChange; 

        for i = 1:length(remainingBranches)
            weight_new(remainingBranches(i)) = weight(remainingBranches(i)) + redistribution(i);
        end
    end

    weight_new = Split(weight, weight_new);

    if ValidateSolution(weight_new, c, V, tolerance) == 0
        weight_new = weight; 
    end
end


function [weight_new] = Split(weight, weight_new)
    if all(weight_new(1:4) < [54.13, 21.56, 34.08, 49.19]) && all(weight_new(1:4) >= 0)
        if weight_new(1) ~= weight(1)
            OK = false;
            while ~OK
                weight_new(5) = rand * weight_new(1);
                weight_new(6) = weight_new(1) - weight_new(5);
                OK = all(weight_new(5:6) < [33.03, 21.84]) && all(weight_new(5:6) >= 0);
            end
        end

        if weight_new(2) ~= weight(2)
            OK = false;
            while ~OK
                weight_new(7) = rand * weight_new(2);
                weight_new(8) = weight_new(2) - weight_new(7);
                OK = all(weight_new(7:8) < [29.96, 24.87]) && all(weight_new(7:8) >= 0);
            end
        end

        if weight_new(4) ~= weight(4)
            OK = false;
            while ~OK
                weight_new(10) = rand * weight_new(4);
                weight_new(9) = weight_new(4) - weight_new(10);
                OK = all(weight_new(9:10) < [47.24, 33.97]) && all(weight_new(9:10) >= 0);
            end
        end

        if any(weight_new([9, 3, 8]) ~= weight([9, 3, 8]))
            OK = false;
            while ~OK
                combinedFlow = weight_new(9) + weight_new(3) + weight_new(8);
                weight_new(11) = rand * combinedFlow;
                weight_new(12) = rand * combinedFlow;
                weight_new(13) = combinedFlow - weight_new(11) - weight_new(12);
                weight_new(17) = weight_new(10) + weight_new(11);

                OK = all(weight_new([11, 12, 13, 17]) < [26.89, 32.76, 39.98, 59.73]) && ...
                          all(weight_new([11, 12, 13, 17]) > 0);
            end
        end

        if any(weight_new([7, 6, 13]) ~= weight([7, 6, 13]))
            OK = false;
            while ~OK
                combinedFlow2 = weight_new(6) + weight_new(7) + weight_new(13);
                weight_new(14) = rand * 37.12;
                weight_new(15) = combinedFlow2 - weight_new(14);
                weight_new(16) = weight_new(5) + weight_new(14);

                OK = all(weight_new([14, 15, 16]) < [37.12, 53.83, 61.65]) && ...
                          all(weight_new([14, 15, 16]) >= 0);
            end
        end

        if any(weight_new([10, 11]) ~= weight([10, 11]))
            weight_new(17) = weight_new(10) + weight_new(11);
        end
    end
end


function [OK] = ValidateSolution(weight, c, V, tolerance)
    OK = false;
    if any(weight > c) || any(weight < 0)
        return; 
    end
    if abs(weight(1) + weight(2) + weight(3) + weight(4) - V) < tolerance && ...
       abs(weight(1) - (weight(5) + weight(6))) < tolerance && ...
       abs(weight(2) - (weight(7) + weight(8))) < tolerance && ...
       abs(weight(4) - (weight(9) + weight(10))) < tolerance && ...
       abs(weight(3) + weight(8) + weight(9) - (weight(11) + weight(12) + weight(13))) < tolerance && ...
       abs(weight(13) + weight(7) + weight(6) - (weight(14) + weight(15))) < tolerance && ...
       abs(weight(14) + weight(5) - weight(16)) < tolerance && ...
       abs(weight(11) + weight(10) - weight(17)) < tolerance && ...
       abs(weight(17) + weight(12) + weight(15) + weight(16) - V) < tolerance
       OK = true; 
    end
end


function [offspring1, offspring2] = PerformCrossover(parent1, parent2, c, V, tolerance)
    offspring1 = zeros(size(parent1));
    offspring2 = zeros(size(parent2));
    offspring1 = (parent1 + parent2) / 2;
    offspring2 = offspring1; 
    if ~ValidateSolution(offspring1, c, V, tolerance) || ~ValidateSolution(offspring2, c, V, tolerance)
        offspring1 = parent1;
        offspring2 = parent2;
    end
end


function [] = RunGeneticAlgorithm()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % If needed remove implicitly set Random Number Generator although it %
    % was left for debugging purposes and standard plot comparisons.      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rng(47);
    a = [1.25 1.25 1.25 1.25 1.25 1.5 1.5 1.5 1.5 1.5 1 1 1 1 1 1 1];
    c = [54.13 21.56 34.08 49.19 33.03 21.84 29.96 24.87 47.24 33.97 ...
         26.89 32.76 39.98 37.12 53.83 61.65 59.73];
    t = ones([1, 17]) * 3;
    V = 100;
    tolerance = 1e-6;

    initialPopulation = zeros(10, 17);
    for i = 1:10
        initialPopulation(i, :) = GenerateInitialPopulation(V);
    end

    fitnessScores = FitnessFunc(a, c, initialPopulation, t, tolerance);
    [selectionCounts, selectedIndices] = PerformRouletteSelection(fitnessScores);
    generationCount = 0;
    fitnessEvolution = zeros(1000, 10);
    fitnessEvolution(1, :) = fitnessScores';

    while true
        generationCount = generationCount + 1;
        newPopulation = GenerateNextPopulation(initialPopulation, selectionCounts, selectedIndices);
        fitnessEvolution(generationCount, :) = FitnessFunc(a, c, newPopulation, t, tolerance)';
        mutationIndex = randi([1, 10]);
        newPopulation(mutationIndex, :) = ApplyMutation(newPopulation(mutationIndex, :), c, V, tolerance);
        [parent1Index, parent2Index] = deal(randi([1, 10]), randi([1, 10]));
        [newPopulation(parent1Index, :), newPopulation(parent2Index, :)] = ...
            PerformCrossover(newPopulation(parent1Index, :), newPopulation(parent2Index, :), c, V, tolerance);
        if generationCount > 1
            change = norm((fitnessEvolution(generationCount - 1, :) - fitnessEvolution(generationCount, :)) ./ fitnessEvolution(generationCount - 1, :));
            if change < 1e-20
                break;
            end
        end
        [selectionCounts, selectedIndices] = PerformRouletteSelection(fitnessEvolution(generationCount, :));
        initialPopulation = newPopulation;
    end

    fprintf("Number of generations until convergence is %d\n", generationCount);
    fprintf("Final Chromosome (NEW POPULATION) is:\n");
    fprintf("%.2f ", newPopulation(10, :));
    fprintf("\n");
    totalTime = 0;
    for j = 1:17
        totalTime = totalTime + t(j) + a(j) * newPopulation(10, j) / (1 - (newPopulation(10, j) / c(j)));
    end
    fprintf("Total time is %.2f\n", totalTime);
    if ValidateSolution(newPopulation(10, :), c, V, tolerance) == true
        fprintf("PASSED.\n");
    else
        fprintf("FAILED\n");
    end
    figure;
    plot(1:generationCount, fitnessEvolution(1:generationCount), '-o');
    title("Fitness Evolution Across Generations STANDARD V");
    xlabel("Generations");
    ylabel("Fitness (Inverse of Total Time)");
    legend("Fitness", 'Location', 'best');
end

RunGeneticAlgorithm();
