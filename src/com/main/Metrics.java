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


    //do Hypervolume
    //reference point for HV - individual with the objectives with the worst values found in set
    //ALGOTYRM Z JMETAL NIE UŻYWA PUNKTU REFERENCYJNEGO, LICZY WZGLĘDEM POCZĄTKU UKŁADU WSPÓŁRZEDNYCH
    private Individual generateReferencePoint(List<Individual> set1, List<Individual> set2, int numOfObjectives,
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


    /*
   returns true if 'point1' dominates 'points2' with respect to the
   to the first 'noObjectives' objectives
   */
    private static boolean dominates(Individual point1, Individual point2, int noOfObjectives) {
        int i;
        int betterInAnyObjective;

        betterInAnyObjective = 0;
        for (i = 0; i < noOfObjectives &&
                point1.getFitnessOfObjectiveByIndex(i) >= point2.getFitnessOfObjectiveByIndex(i); i++) {
            if (point1.getFitnessOfObjectiveByIndex(i) > point2.getFitnessOfObjectiveByIndex(i)) {
                betterInAnyObjective = 1;
            }
        }

        return ((i >= noOfObjectives) && (betterInAnyObjective > 0));
    }

    /*
    Swaps front[i][] with front[j][]
     */
    private static void swap(List<Individual> front, int i, int j) { /* chyba nie używamy, bo zmieniłem funkcję "filterNondominatedSet"
    tak, że nie swapuje, tylko usuwa zdiminowane osobniki */
        Individual temp;

        temp = front.get(i);
        front.set(i, front.get(j));
        front.set(j, temp);
    }
//        private void swapOryginal(double[][] front, int i, int j) {
//        double[] temp;
//
//        temp = front[i];
//        front[i] = front[j];
//        front[j] = temp;
//    }


    /* all nondominated points regarding the first 'noObjectives' dimensions
    are collected; the points referenced by 'front[0..noPoints-1]' are
    considered; 'front' is resorted, such that 'front[0..n-1]' contains
    the nondominated points; n is returned */
    private static int filterNondominatedSet(List<Individual> set, int noPoints, int noOfObjectives) {
        int i = 0, j;
        int n = noPoints;
//        List<Individual> filteredFront = set;

        while (i < n) {
            j = i + 1;
            while (j < n) {
                if (dominates(set.get(i), set.get(j), noOfObjectives)) {
                    /* remove point 'j' */
                    n--;
                    swap(set, j, n);
                } else if (dominates(set.get(j), set.get(i), noOfObjectives)) {
	                /* remove point 'i' */
                    n--;
                    swap(set, i, n);
                    i--;
                    break;
                } else {
                    j++;
                }
            }
            i++;
        }

        return n;
    }
//    private int filterNondominatedSetOryginal(double[][] front, int noPoints, int noObjectives) {
//        int i, j;
//        int n;
//
//        n = noPoints;
//        i = 0;
//        while (i < n) {
//            j = i + 1;
//            while (j < n) {
//                if (dominates(front[i], front[j], noObjectives)) {
//                    /* remove point 'j' */
//                    n--;
//                    swap(front, j, n);
//                } else if (dominates(front[j], front[i], noObjectives)) {
//	/* remove point 'i'; ensure that the point copied to index 'i'
//	   is considered in the next outer loop (thus, decrement i) */
//                    n--;
//                    swap(front, i, n);
//                    i--;
//                    break;
//                } else {
//                    j++;
//                }
//            }
//            i++;
//        }
//        return n;
//    }


    /* calculate next value regarding dimension 'objective'; consider
     points referenced in 'front[0..noPoints-1]' */
    private static double surfaceUnchangedTo(List<Individual> set, int objective) {
        double value;

        if(set.size() >= 1){
            double minValue = set.get(0).getFitnessOfObjectiveByIndex(objective);

            for (Individual point : set) {
                value = point.getFitnessOfObjectiveByIndex(objective);
                if (value < minValue) {
                    minValue = value;
                }
            }

            return minValue;
        }

        //jakiś error czy exception zamist return 0?
        return 0;
    }
//    private double surfaceUnchangedToOryginal(double[][] front, int noPoints, int objective) {
//        int i;
//        double minValue, value;
//
//        if (noPoints < 1) {
//            new JMetalException("run-time error");
//        }
//
//        minValue = front[0][objective];
//        for (i = 1; i < noPoints; i++) {
//            value = front[i][objective];
//            if (value < minValue) {
//                minValue = value;
//            }
//        }
//        return minValue;
//    }

    /* remove all points which have a value <= 'threshold' regarding the
     dimension 'objective'; the points referenced by
     'front[0..noPoints-1]' are considered; 'front' is resorted, such that
     'front[0..n-1]' contains the remaining points; 'n' is returned */
    // trzeba zwrócić uwagę na kryteria MIN, nie wiem czy tutaj je inaczej obsłużyć, czy wcześniej zamienić je na MAX?
    private static int reduceNondominatedSet(List<Individual> set, int noPoints, int objective, double threshold) {
        int n = noPoints;

        for (int i = 0; i < n; i++){
            if(set.get(i).getFitnessOfObjectiveByIndex(objective) <= threshold){
                n--;
                swap(set, i, n);
            }
        }

        return n;
    }
//    private int reduceNondominatedSetOryginal(double[][] front, int noPoints, int objective,
//                                      double threshold) {
//        int n;
//        int i;
//
//        n = noPoints;
//        for (i = 0; i < n; i++) {
//            if (front[i][objective] <= threshold) {
//                n--;
//                swap(front, i, n);
//            }
//        }
//
//        return n;
//    }

    //Hypervolume - HV: ACCURACY, DIVERSITY, CARDINALITY
    public static double HV_metric(List<Individual>  front, int noPoints, int noObjectives) {
        int n = noPoints;
        double volume = 0, distance = 0;

        while(n > 0){
            int nonDominatedPoints;
            double tempVolume, tempDistance;

            nonDominatedPoints = filterNondominatedSet(front, n, noObjectives - 1);

            if (noObjectives < 3) {
                if (nonDominatedPoints < 1) {
                    return -1;
                }
                tempVolume = front.get(0).getFitnessOfObjectiveByIndex(0);
            } else {
                tempVolume = HV_metric(front, n, noObjectives - 1);
            }

            tempDistance = surfaceUnchangedTo(front, noObjectives - 1);
            volume += tempVolume * (tempDistance - distance);
            distance = tempDistance;
            n = reduceNondominatedSet(front, n, noObjectives - 1, distance);
        }

        return volume;
    }
//    public double calculateHypervolumeOryginal(double[][] front, int noPoints, int noObjectives) {
//        int n;
//        double volume, distance;
//
//        volume = 0;
//        distance = 0;
//        n = noPoints;
//        while (n > 0) {
//            int nonDominatedPoints;
//            double tempVolume, tempDistance;
//
//            nonDominatedPoints = filterNondominatedSet(front, n, noObjectives - 1);
//            if (noObjectives < 3) {
//                if (nonDominatedPoints < 1) {
//                    new JMetalException("run-time error");
//                }
//
//                tempVolume = front[0][0];
//            } else {
//                tempVolume = calculateHypervolumeOryginal(front, nonDominatedPoints, noObjectives - 1);
//            }
//
//            tempDistance = surfaceUnchangedTo(front, n, noObjectives - 1);
//            volume += tempVolume * (tempDistance - distance);
//            distance = tempDistance;
//            n = reduceNondominatedSet(front, n, noObjectives - 1, distance);
//        }
//        return volume;
//    }



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
