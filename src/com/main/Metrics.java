package com.main;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class Metrics {

//Metrics reffering to CARDINALITY

    /**
     * @param nonDominatedSet - our examined set of individuals
     * @param trueParetoFront - true known front of Pareto optimal individuals
     * @return Error ratio - number of individuals, that are NOT in true Pareto front / size of set
     */
    //Error Ratio (ER): CARDINALITY ; true PF is needed
    public static double errorRatio(List<Individual> nonDominatedSet, List<Individual> trueParetoFront){
        double ER = 0;
        for (Individual ind: nonDominatedSet) {
            for (Individual optInd: trueParetoFront) {
                if (!ind.equals(optInd)) ER+=1;
            }
        }
        ER = ER/nonDominatedSet.size();
        return ER;
    }

    /**
     * @param nonDominatedSet - set of NON DOMINATED individuals
     * @param popSize - size of whole population
     * @return ratio of non dominated individuals in whole population
     */
    //Ratio of non dominated individuals: CARDINALITY
    public static double ratioOfNonDominated(List<Individual> nonDominatedSet, int popSize){
        double RNI = nonDominatedSet.size()/popSize;
        return RNI;
    }

    /**
     * @param set - set of individuals
     * @return cardinality of set
     */
    //Overall Non-dominated Vector Generation: CARDINALITY
    public static int ONVG(List<Individual> set){
        int ONVG = set.size();
        return ONVG;
    }


//Metrics reffering to ACCURACY
    /**
     * Generetes 1 artifical (artifically optimal) individual with best value on each objective (better than any existing
     * individual).
     * @param set1 - set of solutions from one algorithm
     * @param set2 - set of solutions from another algorithm
     * @param numOfObjectives - number of objectives in mulit-objective EA
     * @param idsOfMinObjective - indexes of objectives that are minimalized (in our problem there is sometimes 1 objective
     *                         minimalized, all others are max); if no minimalized objective - empty list
     * @return Individual optimalIndividual - artificially optimal Individual
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
            for (int i=0; i < numOfObjectives; i++) { //iterating through fitness of objectives of Individual ind
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
            for (int i=0; i < numOfObjectives; i++) { //iterating through fitness of objectives of Individual ind
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

    /**
     * Euclidean distance between two individuals
     * @param a - individual 1
     * @param b - individual 2
     */
    public static double distance(Individual a, Individual b){
        if(a == null || b == null) return -1;
        if(a.getFitnessOfObjectives().size() != b.getFitnessOfObjectives().size()){
//            throw new Exception("Individuals have different number of objectives!");
            System.out.println("Individuals have different number of objectives!");
            return -1;
        }
        if(a.equals(b) || b.equals(a)) return 0;

        double d2 = 0, d;
        for (int i = 0; i < a.getFitnessOfObjectives().size(); i++){
            d2 += Math.pow(a.getFitnessOfObjectiveByIndex(i) - b.getFitnessOfObjectiveByIndex(i), 2);
        }
        d = Math.sqrt(d2);

        return d;
    }

    /**
     *
     * @param set - examined set of individuals
     * @param paretoOptimalSet - true Pareto optimal set or artificial optimal set of individuals (or set of 1 individual)
     *      see Metrics.generateOptimalIndividual()
     *
     * @return Generational distance
     */
    //Generational distance - GD: ACCURACY
    public static double generationalDistance(List<Individual> set, List<Individual> paretoOptimalSet){
        double d = 0, min_d = Double.MAX_VALUE, sumd2 = 0, GD = -1;
        for (Individual ind: set) {//iterate through examined individuals
            for (Individual optimalInd: paretoOptimalSet) {//iterate thorugh optimal individuals, looking for min d^2
                d = distance(ind, optimalInd);
                if (d < min_d) min_d = d;
            }
            sumd2 += Math.pow(min_d, 2);
        }
        GD = Math.sqrt(sumd2)/set.size();

        return GD;
    }

    //stare, był jeden błąd i szukame min_d^2, a nie min_d, chociaż powinno wyjść to samo
//    public static double generationalDistance(List<Individual> set, List<Individual> paretoOptimalSet){
//        double d2 = 0, min_d2 = Double.MAX_VALUE, sumd2 = 0, GD = -1;
//        for (Individual ind: set) {//iterate through examined individuals
//            for (Individual optimalInd: paretoOptimalSet) {//iterate thorugh optimal individuals, looking for min d^2
//                d2 = Math.pow(distance(ind, optimalInd), 2);
//                if (d2 < min_d2) min_d2 = d2;
//            }
//            sumd2 += min_d2; // += min_d2 !!!
//        }
//        GD = Math.sqrt(sumd2)/set.size();
//
//        return GD;
//    }

    /**
     *
     * @param set - examined set of individuals
     * @param paretoOptimalSet - true Pareto optimal set or artificall optimal set of individuals (or set of 1 individual)
     *      see Metrics.generateOptimalIndividual()
     *
     * @return Inverted Generational distance
     */
    //Inverted Generational distance - IGD: ACCURACY (also diversity, but not in our case)
    public static double invertedGenerationalDistance(List<Individual> set, List<Individual> paretoOptimalSet){
        double IGD = -1, d = 0, sum_d = 0, min_d = Double.MAX_VALUE;

        for (Individual ind: set) {//iterate thorugh test individuals
            for (Individual optimalInd: paretoOptimalSet) {//iterate through optimal individuals, looking for min d^2
                d = distance(ind, optimalInd);
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
        double MPFE = 0, d = 0, min_d = Double.MAX_VALUE, max_d = -1;
        List<Double> distancesToNearestOptimal = new ArrayList<>();

        for (Individual ind: set) {//iterate thorugh test individuals
            d = 0;
            for (Individual optimalInd: paretoOptimalSet) {//iterate thorugh optimal individuals, looking for min d for this ind
//                for (int i = 0; i < ind.getFitnessOfObjectives().size(); i++) {//calculate d^2
//                    d += Math.pow(ind.getFitnessOfObjectiveByIndex(i) - optimalInd.getFitnessOfObjectiveByIndex(i), 2);
//                }
//                d = Math.sqrt(d);
//                if (d > max_d) max_d = d;
                d = distance(ind, optimalInd);
                if (d < min_d) min_d = d; // odległość do najbliższego pkt z paretoOptimalSet (min)
            }
            distancesToNearestOptimal.add(min_d);
        }
        Collections.sort(distancesToNearestOptimal);// sorting ascending
        MPFE = distancesToNearestOptimal.get(distancesToNearestOptimal.size() - 1);// getting last (max) number

        return MPFE;
    }


    //funkcje do Hypervolume:

    //reference point for HV - individual with the objectives with the worst values found in set
    //ALGOTYRM Z JMETAL NIE UŻYWA PUNKTU REFERENCYJNEGO, LICZY WZGLĘDEM POCZĄTKU UKŁADU WSPÓŁRZĘDNYCH
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

        // TODO jak byśmy używali, to trzeba te wartości jeszcze "pogorszyć", żeby nie były równe najgorszym
        // TODO ze znalezionych, ale jeszcze gorsze od nich
        for (int i = 0; i < worstArray.length; i++) {
            referencePoint.getFitnessOfObjectives().add(worstArray[i]);
        }

        return referencePoint;
    }

    /*
   returns true if 'point1' dominates 'points2' with respect to the
   to the first 'noObjectives' objectives
   */
    /*private*/public static boolean dominates(Individual point1, Individual point2, int noOfObjectives) { // działa
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
//    private boolean dominatesOryginal(double point1[], double point2[], int noObjectives) {
//        int i;
//        int betterInAnyObjective;
//
//        betterInAnyObjective = 0;
//        for (i = 0; i < noObjectives && point1[i] >= point2[i]; i++) {
//            if (point1[i] > point2[i]) {
//                betterInAnyObjective = 1;
//            }
//        }
//
//        return ((i >= noObjectives) && (betterInAnyObjective > 0));
//    }

    /*
    Swaps individual [i] with individual [j] in list "front"
     */
    /*private*/public static void swap(List<Individual> front, int i, int j) { // działa
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


    /* all nondominated points regarding the first 'noObjectives' dimensions are collected; the points referenced by
    'front[0..noPoints-1]' are considered; 'front' is resorted, such that 'front[0..n-1]' contains
    the nondominated points; n is returned */
    /*private*/public static int filterNondominatedSet(List<Individual> set, int noPoints, int noOfObjectives) { // działa
        int i = 0, j;
        int n = noPoints;

        while (i < n) {
            j = i + 1;
            while (j < n) {
                if (dominates(set.get(i), set.get(j), noOfObjectives)) {
                    /* remove point 'j' */
                    n--;
                    swap(set, j, n);
                } else if (dominates(set.get(j), set.get(i), noOfObjectives)) {
                /* remove point 'i'; ensure that the point copied to index 'i'
                   is considered in the next outer loop (thus, decrement i) */
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


    /* calculate min value regarding dimension 'objective'; consider
     points referenced in 'front[0..noPoints-1]'
     znajduje minimalną wartość na danym kryterium spośród wszytstkich osobników*/
    /*private*/public static double surfaceUnchangedTo(List<Individual> set, int objective) { // działa
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

        //jakiś error czy exception zamist return -1?
        return -1;
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

    /* remove all points which have a value <= 'threshold' regarding the dimension 'objective'; the points referenced by
     'front[0..noPoints-1]' are considered; 'front' is resorted, such that 'front[0..n-1]' contains the remaining points;
     'n' is returned */
    // TODO trzeba zwrócić uwagę, gdy kryteria są MIN (a nie MAX), nie wiem czy tutaj je inaczej obsłużyć,
    // TODO czy gdzieś wcześniej zamienić je na MAX?
    /*private*/public static int reduceNondominatedSet(List<Individual> set, int noPoints, int objective, double threshold) {
        int n = noPoints;
        int i;

        for (i = 0; i < n; i++){
            if(set.get(i).getFitnessOfObjectiveByIndex(objective) <= threshold){
                n--;
                swap(set, i, n);
                i--;// DODANA LINIA !!!
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
    public static double calculateHypervolume(List<Individual> front, int noPoints, int noObjectives) {
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
                tempVolume = calculateHypervolume(front, nonDominatedPoints, noObjectives - 1);
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

//    /**
//     * @param a Individual a
//     * @param b Individual b
//     * @return euclidean distance between Individuals a and b
//     */
//    public static double distance(Individual a, Individual b){
//        double d2= 0, d = -1;
//        if (a.getFitnessOfObjectives().size() != b.getFitnessOfObjectives().size()) return -1;
//        for (int i = 0; i < a.getFitnessOfObjectives().size(); i++) {//calculate d
//            d2 += Math.pow(a.getFitnessOfObjectiveByIndex(i) - b.getFitnessOfObjectiveByIndex(i), 2);
//        }
//        d = Math.sqrt(d2);
//
//        return d;
//    }

    // funkcje do Uniform Distribution:
    public static int sh(Individual a, Individual b, double sigma){
        if (distance(a,b) < sigma) return 1;
        else return 0;
    }

    public static int nicheCount(Individual a, List<Individual> set, double sigma){
        int count = 0;
        for (Individual ind: set ) {
            if (!(ind == a)){ //if("ind" nie jest tym samym obiektem, co "a") - liczymy sh() z każdym osobnikiem, oprócz siebie samego
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



    // funkcja do Spacing:
    //mean of distances to closest point from true Pareto optimal front
    public static double meanDistance(List<Individual> set, List<Individual> paretoOptimalSet){
        double d = 0, min_d = Double.MAX_VALUE, sum_d = 0, mean_d = 0;

//        for (Individual ind: set) {//iterate thorugh test individuals
//            for (Individual optimalInd: paretoOptimalSet) {//iterate thorugh optimal individuals, summing d
//                d2 = 0;
//                for (int i = 0; i < ind.getFitnessOfObjectives().size(); i++) {//calculate d^2
//                    d2 += Math.pow(ind.getFitnessOfObjectiveByIndex(i) - optimalInd.getFitnessOfObjectiveByIndex(i), 2);
//                }
//                d = Math.sqrt(d2);
//                sum_d += d;
//            }
//            mean_d = sum_d/paretoOptimalSet.size();
//            wholeMean_d += mean_d;
//        }
//        wholeMean_d = wholeMean_d/set.size();
        for (Individual ind: set) {//iterate thorugh test individuals
            min_d = Double.MAX_VALUE;
            for (Individual optimalInd: paretoOptimalSet) {//iterate thorugh optimal individuals, searching for min d
                d = distance(ind, optimalInd);
                if(d < min_d) min_d = d;
            }
            mean_d += min_d;
        }
        mean_d = mean_d/set.size();// mean of  sum(min_d[i]), i = 1,..., size_of_set

        return mean_d;
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

// Number of Distinct Choices (NDC) - można zaimplementować, jak nam się będzie chiciało xd



//SPREAD
    //Maximum Spread (MS)
    public static double maximumSpread(List<Individual> set, List<Individual> paretoOptimalSet, int dimensions){
        //double max_set = -Double.MAX_VALUE,  min_set = Double.MAX_VALUE,  max_optimal = -Double.MAX_VALUE,  min_optimal = Double.MAX_VALUE;
        double[] max_set = new double[dimensions];
        double[] min_set = new double[dimensions];
        double[] max_trueParetoFront = new double[dimensions];
        double[] min_trueParetoFront = new double[dimensions];
        Arrays.fill(max_set, -Double.MAX_VALUE);
        Arrays.fill(max_trueParetoFront, -Double.MAX_VALUE);
        Arrays.fill(min_set, Double.MAX_VALUE);
        Arrays.fill(min_trueParetoFront, Double.MAX_VALUE);

        double[] min = new double[dimensions];
        double[] max = new double[dimensions];
        double sum = 0, MS = -1;

        for (int i = 0; i < dimensions; i++) {
            for (Individual ind: set) {
                if (ind.getFitnessOfObjectiveByIndex(i) > max_set[i]) max_set[i] = ind.getFitnessOfObjectiveByIndex(i);
                if (ind.getFitnessOfObjectiveByIndex(i) < min_set[i]) min_set[i] = ind.getFitnessOfObjectiveByIndex(i);
            }
            for (Individual ind: paretoOptimalSet) {
                if (ind.getFitnessOfObjectiveByIndex(i) > max_trueParetoFront[i]) max_trueParetoFront[i] = ind.getFitnessOfObjectiveByIndex(i);
                if (ind.getFitnessOfObjectiveByIndex(i) < min_trueParetoFront[i]) min_trueParetoFront[i] = ind.getFitnessOfObjectiveByIndex(i);
            }
        }

        for (int i = 0; i < dimensions; i++) {
            max[i] = Math.min(max_set[i],max_trueParetoFront[i]);
            min[i] = Math.max(min_set[i],min_trueParetoFront[i]);

            sum += Math.pow((max[i] - min[i]) / (max_trueParetoFront[i] - min_trueParetoFront[i]),2);
        }
        MS = Math.sqrt(sum/dimensions);

        return MS;
    }

}
