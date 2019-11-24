package com.main;

import java.util.ArrayList;
import java.util.List;

public class Individual {

    private List<Double> fitnessOfObjectives = new ArrayList<>();


    public Individual(List<Double> fitnessOfObjectives) {
        this.fitnessOfObjectives = fitnessOfObjectives;
    }

    public Individual() { }


    public List<Double> getFitnessOfObjectives() {
        return fitnessOfObjectives;
    }

    public void setFitnessOfObjectives(List<Double> fitnessOfObjectives) {
        this.fitnessOfObjectives = fitnessOfObjectives;
    }

    public double getFitnessOfObjectiveByIndex(int id) {
        return fitnessOfObjectives.get(id);
    }

    @Override
    public boolean equals(Object b){
        if (this == b) return true;
        if (b == null) return false;
        if (this.getClass() != b.getClass()) return false;
        Individual other = (Individual) b;
        if (fitnessOfObjectives.size() == other.getFitnessOfObjectives().size()){
            for (int i = 0; i < fitnessOfObjectives.size(); i++) {
                if(fitnessOfObjectives.get(i) != null && other.getFitnessOfObjectives().get(i) != null){
                    if(fitnessOfObjectives.get(i) != other.getFitnessOfObjectiveByIndex(i))
                        return false;
                } else return false;
            }
        }else return false;

        return true;
    }
}
