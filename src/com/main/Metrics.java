package com.main;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class Metrics {

//Metrics reffering to CARDINALITY

    //Error Ratio (ER): CARDINALITY ; true PF is needed
    public static double errorRatio(List<Individual> nonDominatedSet, List<Individual> paretoOptimalSet){
        double ER = 0;
        for (Individual ind: nonDominatedSet) {
            for (Individual optInd: paretoOptimalSet) {
                if (!ind.equals(optInd)) ER+=1;
            }
        }
        ER = ER/nonDominatedSet.size();
        return ER;
    }

    /**
     * @param nonDominatedSet
     * @param popSize - size of whole population
     * @return ratio of non dominated individuals in whole population
     */
    //Ratio of non dominated individuals: CARDINALITY
    public static double ratioOfNonDominated(List<Individual> nonDominatedSet, int popSize){
        double RNI = nonDominatedSet.size()/popSize;
        return RNI;
    }

    //Overall Non-dominated Vector Generation: CARDINALITY
    public static int ONVG(List<Individual> set){
        int ONVG = set.size();
        return ONVG;
    }


//Metrics reffering to ACCURACY
    /**
     * Generetes 1 artifical (artifically optimal) individual with best value on each objective (better than existing
     * individuals).
     * @param set1 - set of solutions from one algorithm
     * @param set2 - set of solutions from another algorithm
     * @param numOfObjectives - number of objectives in mulit-objective EA
     * @param idsOfMinObjective - index of objective that is minimalized (in our problem there is sometimes 1 objective
     *                         minimalized, all others are max); if no minimalized objective - e.g. idOfMinObjective = -1
     * @return Individual optimalIndividual - artifically optimal Individual
     */
    public static Individual generateOptimalIndividual(List<Individual> set1, List<Individual> set2, int numOfObjectives,
                                                List<Integer> idsOfMinObjective){
        Individual optimalIndividual = new Individual();
        double bestArray [] = new double[numOfObjectives];//best value of each objective
        //bestArray[idOfMinObjective] = Double.MAX_VALUE;
        for (int i=0; i<numOfObjectives; i++) //all objectives as max
            bestArray[i] = -1;
        for (int i: idsOfMinObjective) //only objectives from idsOfMinObjective as min
            bestArray[i] = Double.MAX_VALUE;

        //searching for best values of each obiective in set1 and set2
        for (Individual ind: set1) {//iterating through set1 individuals
            for (int i=0; i<ind.getFitnessOfObjectives().size(); i++) { //iterating through fitness of objectives of Individual ind
                if (idsOfMinObjective.contains(i)){//min
                    if (ind.getFitnessOfObjectiveByIndex(i) < bestArray[i])
                        bestArray[i] = ind.getFitnessOfObjectiveByIndex(i);
                } else{//max
                    if (ind.getFitnessOfObjectiveByIndex(i) > bestArray[i])
                        bestArray[i] = ind.getFitnessOfObjectiveByIndex(i);
                }
            }
        }
        for (Individual ind: set2) {//iterating through set2 individuals
            for (int i=0; i<ind.getFitnessOfObjectives().size(); i++) { //iterating through fitness of objectives of Individual ind
                if (idsOfMinObjective.contains(i)){//min
                    if (ind.getFitnessOfObjectiveByIndex(i) < bestArray[i])
                        bestArray[i] = ind.getFitnessOfObjectiveByIndex(i);
                } else{//max
                    if (ind.getFitnessOfObjectiveByIndex(i) > bestArray[i])
                        bestArray[i] = ind.getFitnessOfObjectiveByIndex(i);
                }
            }
        }

        //upgrading each objective for optimalIndividual by 100%
        for (int i = 0; i < bestArray.length; i++) {
            if (idsOfMinObjective.contains(i)) //min
                bestArray[i] =  bestArray[i]/2;
            else //max
                bestArray[i] =  bestArray[i]*2;
        }

        for (int i = 0; i < bestArray.length; i++) {
            optimalIndividual.getFitnessOfObjectives().add(bestArray[i]);
        }

        return optimalIndividual;
    }


    //Generational distance - GD: ACCURACY
    public static double generationalDistance(List<Individual> set, List<Individual> paretoOptimalSet){
        double d2 = 0, min_d2 = Double.MAX_VALUE, sumd2 = 0, GD = 0;
        for (Individual ind: set) {//iterate thorugh test individuals
            for (Individual optimalInd: paretoOptimalSet) {//iterate thorugh optimal individuals, looking for min d^2
                d2 = 0;
                for (int i = 0; i < ind.getFitnessOfObjectives().size(); i++) {//calculate d^2
                    d2 += Math.pow(ind.getFitnessOfObjectiveByIndex(i) - optimalInd.getFitnessOfObjectiveByIndex(i), 2);
                }
                if (d2 < min_d2) min_d2 = d2;
            }
            sumd2 += d2;
        }
        GD = Math.sqrt(sumd2)/set.size();

        return GD;
    }


    //Inverted Generational distance - IGD: ACCURACY (also diversity, but not in our case)
    public static double invertedGenerationalDistance(List<Individual> set, List<Individual> paretoOptimalSet){
        double IGD = -1, d2 = 0, d = 0, sum_d = 0, min_d = Double.MAX_VALUE;
        /*
        for (Individual ind: set) {
            d = 0;
            for (int i = 0; i < ind.getFitnessOfObjectives().size(); i++) {//calculate d between optimalIndividual and individual from set
                d += Math.pow(ind.getFitnessOfObjectiveByIndex(i) - optimalIndividual.getFitnessOfObjectiveByIndex(i), 2);
            }
            d = Math.sqrt(d);
            //check if actual d is a new min
            if (d < min_d) min_d = d;
        }
        */
        for (Individual ind: set) {//iterate thorugh test individuals
            for (Individual optimalInd: paretoOptimalSet) {//iterate thorugh optimal individuals, looking for min d^2
                d2 = 0;
                for (int i = 0; i < ind.getFitnessOfObjectives().size(); i++) {//calculate d^2
                    d2 += Math.pow(ind.getFitnessOfObjectiveByIndex(i) - optimalInd.getFitnessOfObjectiveByIndex(i), 2);
                }
                d = Math.sqrt(d2);
                if (d < min_d) min_d = d;
            }
            sum_d += d;
        }
        IGD = sum_d/paretoOptimalSet.size();

        return IGD;
    }


    /*Maximum Pareto Front Error - MPFE: ACCURACY
    It measures the largest distance in the objective space between any individual in the
    approximation front and the corresponding closest vector in the true Pareto front.*/
    public static double maximumParetoFrontError(List<Individual> set, List<Individual> paretoOptimalSet){
        double MPFE = 0, d = 0, max_d = -1;
        for (Individual ind: set) {//iterate thorugh test individuals
            for (Individual optimalInd: paretoOptimalSet) {//iterate thorugh optimal individuals, looking for min d^2
                d = 0;
                for (int i = 0; i < ind.getFitnessOfObjectives().size(); i++) {//calculate d^2
                    d += Math.pow(ind.getFitnessOfObjectiveByIndex(i) - optimalInd.getFitnessOfObjectiveByIndex(i), 2);
                }
                d = Math.sqrt(d);
                if (d > max_d) max_d = d;
            }
        }
        MPFE = max_d;

        return MPFE;
    }



    //?not completed? reference point for HV - individual with the objectives with the worst values found in set
    public static Individual generateReferencePoint(List<Individual> set1, List<Individual> set2, int numOfObjectives,
                                                       List<Integer> idsOfMinObjective){
        Individual referencePoint = new Individual();
        double worstArray [] = new double[numOfObjectives];//best value of each objective
        for (int i=0; i<numOfObjectives; i++) //all objectives as max
            worstArray[i] = Double.MAX_VALUE;
        for (int i: idsOfMinObjective) //only objectives from idsOfMinObjective as min
            worstArray[i] = -1;

        //searching for WORST values of each obiective in set1 and set2
        for (Individual ind: set1) {//iterating through set1 individuals
            for (int i=0; i<ind.getFitnessOfObjectives().size(); i++) { //iterating through fitness of objectives of Individual ind
                if (idsOfMinObjective.contains(i)){//looking for the worst value on MIN objective
                    if (ind.getFitnessOfObjectiveByIndex(i) > worstArray[i])
                        worstArray[i] = ind.getFitnessOfObjectiveByIndex(i);
                } else{//looking for the worst value on MAX objectives
                    if (ind.getFitnessOfObjectiveByIndex(i) < worstArray[i])
                        worstArray[i] = ind.getFitnessOfObjectiveByIndex(i);
                }
            }
        }
        for (Individual ind: set2) {//iterating through set2 individuals
            for (int i=0; i<ind.getFitnessOfObjectives().size(); i++) { //iterating through fitness of objectives of Individual ind
                if (idsOfMinObjective.contains(i)){//looking for the worst value on MIN objective
                    if (ind.getFitnessOfObjectiveByIndex(i) > worstArray[i])
                        worstArray[i] = ind.getFitnessOfObjectiveByIndex(i);
                } else{//looking for the worst value on MAX objectives
                    if (ind.getFitnessOfObjectiveByIndex(i) < worstArray[i])
                        worstArray[i] = ind.getFitnessOfObjectiveByIndex(i);
                }
            }
        }

        for (int i = 0; i < worstArray.length; i++) {
            referencePoint.getFitnessOfObjectives().add(worstArray[i]);
        }

        return referencePoint;
    }


    //Hypervolume - HV: ACCURACY, DIVERSITY, CARDINALITY
    public static double HV_metric(List<Individual> set, Individual referencePoint){
        double HV = -1;
        List<Individual> countedIndividuals = new ArrayList<>();
        countedIndividuals.add(referencePoint);

        for (Individual ind: set) {

        }

        return ;
    }


//Metrics reffering to DISTRIBUTION

    /**
     * @param a Individual a
     * @param b Individual b
     * @return euclidean distance between Individuals a and b
     */
    public static double distance(Individual a, Individual b){
        double d2= 0, d = -1;
        if (a.getFitnessOfObjectives().size() != b.getFitnessOfObjectives().size()) return -1;
        for (int i = 0; i < a.getFitnessOfObjectives().size(); i++) {//calculate d
            d2 += Math.pow(a.getFitnessOfObjectiveByIndex(i) - b.getFitnessOfObjectiveByIndex(i), 2);
        }
        d = Math.sqrt(d2);

        return d;
    }

    public static int sh(Individual a, Individual b, double sigma){
        if (distance(a,b) < sigma) return 1;
        else return 0;
    }

    public static int nicheCount(Individual a, List<Individual> set, double sigma){
        int count = 0;
        for (Individual ind: set ) {
            if (!ind.equals(a)){
                count += sh(a, ind, sigma);
            }
        }
        return count;
    }

    public static double meanNicheCount(List<Individual> set, double sigma){
        double sum = 0;
        for (Individual ind: set) {
            sum += nicheCount(ind, set, sigma);
        }
        double mean = sum/set.size();

        return mean;
    }

    public static double standardDeviationOfNicheCount(List<Individual> set, double sigma){
        double std = 0;
        double mean = meanNicheCount(set,sigma);
        for (Individual ind: set) {
            std += Math.pow(nicheCount(ind, set, sigma) - mean, 2);
        }
        std = Math.sqrt(std/(set.size()-1));

        return std;
    }

    //Uniform Distribution (UD)
    public static double uniformDistribution(List<Individual> set, double sigma){
        double UD = -1;
        double std = standardDeviationOfNicheCount(set, sigma);
        UD = 1/(1+std);

        return UD;
    }



    //mean distance of set of individuals to pareto front
    public static double meanDistance(List<Individual> set, List<Individual> paretoOptimalSet){
        double d2 = 0, d = 0, sum_d = 0, mean_d = 0, wholeMean_d = 0;

        for (Individual ind: set) {//iterate thorugh test individuals
            for (Individual optimalInd: paretoOptimalSet) {//iterate thorugh optimal individuals, summing d
                d2 = 0;
                for (int i = 0; i < ind.getFitnessOfObjectives().size(); i++) {//calculate d^2
                    d2 += Math.pow(ind.getFitnessOfObjectiveByIndex(i) - optimalInd.getFitnessOfObjectiveByIndex(i), 2);
                }
                d = Math.sqrt(d2);
                sum_d += d;
            }
            mean_d = sum_d/paretoOptimalSet.size();
            wholeMean_d += mean_d;
        }
        wholeMean_d = wholeMean_d/set.size();

        return wholeMean_d;
    }

    //Spacing
    public static double spacing(List<Individual> set, List<Individual> paretoOptimalSet){
        double d2 = 0, d = 0, min_d2 = Double.MAX_VALUE, sum_d = 0, mean_d = 0, s = -1;
        mean_d = meanDistance(set, paretoOptimalSet);
        for (Individual ind: set) {//iterate thorugh test individuals
            for (Individual optimalInd: paretoOptimalSet) {//iterate thorugh optimal individuals, looking for min d^2
                d2 = 0;
                for (int i = 0; i < ind.getFitnessOfObjectives().size(); i++) {//calculate d^2
                    d2 += Math.pow(ind.getFitnessOfObjectiveByIndex(i) - optimalInd.getFitnessOfObjectiveByIndex(i), 2);
                }
                if (d2 < min_d2) min_d2 = d2;
            }
            d = Math.sqrt(d2);
            sum_d += Math.pow(d - mean_d, 2);
        }
        s = Math.sqrt(sum_d/set.size());

        return s;
    }

// Number of Distinct Choices (NDC)


//SPREAD
    //Maximum Spread (MS)
    public static double maximumSpread(List<Individual> set, List<Individual> paretoOptimalSet, int dimensions){
        //double max_set = -Double.MAX_VALUE,  min_set = Double.MAX_VALUE,  max_optimal = -Double.MAX_VALUE,  min_optimal = Double.MAX_VALUE;
        double[] max_set = new double[dimensions];
        double[] min_set = new double[dimensions];
        double[] max_optimal = new double[dimensions];
        double[] min_optimal = new double[dimensions];
        Arrays.fill(max_set, -Double.MAX_VALUE);
        Arrays.fill(max_optimal, -Double.MAX_VALUE);
        Arrays.fill(min_set, Double.MAX_VALUE);
        Arrays.fill(min_optimal, Double.MAX_VALUE);

        double[] min = new double[dimensions];
        double[] max = new double[dimensions];
        double sum = 0, MS = -1;

        for (int i = 0; i < dimensions; i++) {
            for (Individual ind: set) {
                if (ind.getFitnessOfObjectiveByIndex(i) > max_set[i]) max_set[i] = ind.getFitnessOfObjectiveByIndex(i);
                if (ind.getFitnessOfObjectiveByIndex(i) < min_set[i]) min_set[i] = ind.getFitnessOfObjectiveByIndex(i);
            }
            for (Individual ind: paretoOptimalSet) {
                if (ind.getFitnessOfObjectiveByIndex(i) > max_optimal[i]) max_optimal[i] = ind.getFitnessOfObjectiveByIndex(i);
                if (ind.getFitnessOfObjectiveByIndex(i) < min_optimal[i]) min_optimal[i] = ind.getFitnessOfObjectiveByIndex(i);
            }
        }

        for (int i = 0; i < dimensions; i++) {
            max[i] = Math.min(max_set[i],max_optimal[i]);
            min[i] = Math.max(min_set[i],min_optimal[i]);

            sum += Math.pow((max[i] - min[i]) / (max_optimal[i] - min_optimal[i]),2);
        }
        MS = Math.sqrt(sum/dimensions);

        return MS;
    }

}
