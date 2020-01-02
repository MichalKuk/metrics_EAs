package com.main;

import javax.lang.model.type.NullType;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;


public class Main {

    public static void main(String[] args) {
        List<Individual> individuals = new ArrayList<>();
        List<Individual> trueParetoFront = new ArrayList<>();

        //TESTY
        int noOfObjectives = 3;
        List<Double> l = new ArrayList<>(List.of(10.0, 5.0, 3.5));
        Individual ind = new Individual(l);
        individuals.add(ind);

//        ind = new Individual(l); //na wszelki wypadek, żeby to był inny obiekt o takich samych wartościach
//        trueParetoFront.add(ind);

        l = new ArrayList<>(List.of(5.5, 3.0, 2.0));
        ind = new Individual(l);
        individuals.add(ind);

        l = new ArrayList<>(List.of(3.0, 10.0, 5.5));
        ind = new Individual(l);
        individuals.add(ind);

        l = new ArrayList<>(List.of(3.0, 10.0, 6.0));
        ind = new Individual(l);
        individuals.add(ind);

        l = new ArrayList<>(List.of(10.0, 10.5, 6.0));
        ind = new Individual(l); //na wszelki wypadek, żeby to był inny obiekt o takich samych wartościach
        trueParetoFront.add(ind);

        l = new ArrayList<>(List.of(8.0, 8.0, 9.0));
        ind = new Individual(l); //na wszelki wypadek, żeby to był inny obiekt o takich samych wartościach
        trueParetoFront.add(ind);

//        System.out.println("Individual toString: " + ind);

//        //CARDINALITY
//        System.out.println("Error ratio: " + Metrics.errorRatio(individuals, trueParetoFront));
//        System.out.println("Ratio of non dominated individuals: " + Metrics.ratioOfNonDominated(individuals, 3));
//        System.out.println("ONVG: " + Metrics.ONVG(individuals));
//
//        //ACCURACY
//        Individual OptimalInd = Metrics.generateOptimalIndividual(individuals, Collections.emptyList(),  3, Collections.emptyList());
//        trueParetoFront = new ArrayList<>(List.of(OptimalInd));
//        System.out.println("Optimal individual: " + OptimalInd);
//        System.out.println("Euclidean distance between " + individuals.get(0) + " and " + individuals.get(1) + " : "
//                + Metrics.distance(individuals.get(0), individuals.get(1)));
//
//        System.out.println("Generational Distance: " + Metrics.generationalDistance(individuals, trueParetoFront));
//        System.out.println("Inverted Generational Distance: " + Metrics.invertedGenerationalDistance(individuals, trueParetoFront));
//        System.out.println("Maximum Pareto Front Error: " + Metrics.maximumParetoFrontError(individuals, trueParetoFront));

//        System.out.println("Hypervolume: " + Metrics.calculateHypervolume(individuals, individuals.size(), 3));


        //getFitnessOfObjectiveByIndex(i)
//        System.out.println("Fitness objectives of ind (using getFitnessOfObjectiveByIndex(i)):");
//        for(int i=0; i < ind2.getFitnessOfObjectives().size(); i++){
//            System.out.println(ind2.getFitnessOfObjectiveByIndex(i));
//        }

        //Metrics.dominates()
//        System.out.println("Domination: ");
//        for (Individual i1: individuals) {
//            for (Individual i2: individuals) {
//                if (Metrics.dominates(i1, i2, noOfObjectives)){
//                    System.out.println(i1+" dominates "+i2);
//                } else {
//                    System.out.println(i1+" DOES NOT dominate "+i2);
//                }
//            }
//        }

        //Metrics.swap()
//        int i = 0, j = individuals.size()-1;
//        for(int ii = 0; ii < individuals.size(); ii++){
//            System.out.println(ii+". "+individuals.get(ii));
//        }
//        System.out.println("Swapping " + individuals.get(i) + " with " + individuals.get(j) +":");
//        Metrics.swap(individuals, i, j);
//        for(int ii = 0; ii < individuals.size(); ii++){
//            System.out.println(ii+". "+individuals.get(ii));
//        }

        //Metrics.filterNondominatedSet()
//        System.out.println("\nMetrics.filterNondominatedSet(): ");
//        int nondominatedSetSize;
//        nondominatedSetSize = Metrics.filterNondominatedSet(individuals, individuals.size(), noOfObjectives);
//        for(int ii = 0; ii < individuals.size(); ii++){
//            if(ii < nondominatedSetSize) {
//                System.out.println("NONdominated: " + individuals.get(ii));
//            } else {
//                System.out.println("Dominated: " + individuals.get(ii));
//            }
//        }

        //Metrics.surfaceUnchangedTo()
        for (Individual i1: individuals) {
            System.out.println(i1);
        }
//        for(int ii = 0; ii < noOfObjectives; ii++){
//            System.out.println("Min on objective " + ii + ": " + Metrics.surfaceUnchangedTo(individuals, ii));
//        }

        //Metrics.reduceNondominatedSet()
//        System.out.println("Metrics.reduceNondominatedSet: objective = 0, threshold = 4");
//        int reduced_n = Metrics.reduceNondominatedSet(individuals, individuals.size(), 0, 4.0);
//        System.out.println("reduced_n = " + reduced_n);
//        for(int ii = 0; ii < reduced_n; ii++){
//            System.out.println(individuals.get(ii));
//        }

        //HV
//        double HV = Metrics.calculateHypervolume(individuals, individuals.size(), noOfObjectives);
//        System.out.println("Hypervolume = " + HV);

        //DISTRIBUTION
        //Metrics.uniformDistribution()
//        double UD = Metrics.uniformDistribution(individuals,4.0);
//        System.out.println("Uniform distribution = " + UD);

        //Metrics.spacing()
//        double s = Metrics.spacing(individuals, trueParetoFront);
//        System.out.println("Spacing = " + s);

        //SPREAD
        //Metrics.maximumSpread()
//        double MS = Metrics.maximumSpread(individuals, trueParetoFront, noOfObjectives);
//        System.out.println("Maximum Spread = " + MS);

    }
}
