package crf;

public class VariableNode extends Node {

	
    public int  y;
    public int label_type;
    public double[] state_factor;
    public double[] marginal;
	public VariableNode() {
		// TODO Auto-generated constructor stub
		super();
	}
	
	public void Init(int num_label)
	{
		   BasicInit(num_label);
		   nodetype =0;
		   state_factor = GetDoubleArr(num_label);
		   marginal = GetDoubleArr(num_label);
	}
	public void BeliefPropagation(globalvar gl, boolean labeled_given)
	{
	    double product;
	    for (int i = 0; i < neighbor.size(); i ++)
	    {
	        FactorNode f;
	        f= (FactorNode) neighbor.get(i);
	        for (int y = 0; y < num_label; y ++)
	        {
	            product = state_factor[y];
	            for (int j = 0; j < neighbor.size(); j ++)
	                if (i != j)
	                {
	                	double bpy = ((double[])belief.get(j))[y];
	                    product *= bpy;
	                }
	            msg[y] = product;
	        }

	        NormalizeMessage();
	        f.GetMessageFrom(nodeid, msg, gl);
	    }
	}
	public void MaxSumPropagation(globalvar gl, boolean labeled_given)
	{
	    double product;

	    for (int i = 0; i < neighbor.size(); i ++)
	    {
	        FactorNode f = (FactorNode) neighbor.get(i);
	        for (int y = 0; y < num_label; y ++)        
	        {
	            product = state_factor[y];
	            for (int j = 0; j < neighbor.size(); j ++)
	                if (i != j)
	                {
	                	double bpy = ((double[])belief.get(j))[y];         	
	                    product *= bpy;
	                }
	            msg[y] = product;
	        }
	        NormalizeMessage();
	        f.GetMessageFrom(nodeid, msg, gl);
	    }
	}
}
