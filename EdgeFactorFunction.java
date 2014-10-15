package crf;

import java.util.HashMap;

public class EdgeFactorFunction {

	
	public  int num_label;
	public double[] lambda;
	public HashMap<Integer, Integer> feature_offset;

	public EdgeFactorFunction(int num_label, double[] p_lambda, HashMap<Integer, Integer> edge_feature_offset)
	 {
	        this.num_label = num_label;
	        this.lambda = p_lambda;
	        this.feature_offset = edge_feature_offset;
	}
	    
	 public double GetValue(int y1, int y2)
	  {
	        int a = y1 < y2 ? y1 : y2;
	        int b = y1 > y2 ? y1 : y2;
	        int i = feature_offset.get(a*num_label + b);
	        return Math.exp(lambda[i]);
	  }
}