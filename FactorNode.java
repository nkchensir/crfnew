package crf;

public class FactorNode extends Node {

	public FactorNode() {
		// TODO Auto-generated constructor stub
	}
    public EdgeFactorFunction  func;
    public double[][] marginal;

    public void Init(int num_label)
    {
        BasicInit(num_label);
        nodetype = 1;
        marginal = new double[num_label][];
        for (int i = 0; i < num_label; i ++)
            marginal[i] = new double[num_label];
    }
    
    public void BeliefPropagation(globalvar gl, boolean labeled_given)
    {
        for (int i = 0; i < 2; i ++)
        {
            if (labeled_given && ((VariableNode) neighbor.get(i)).label_type == 0)
            {
                for (int y = 0; y < num_label; y ++)
                    msg[y] = 0;
                msg[((VariableNode) neighbor.get(i)).y] = 1.0;
            }
            else
            {
                for (int y = 0; y < num_label; y ++)
                {
                    double s = 0;
                    for (int y1 = 0; y1 < num_label; y1 ++)
                    {
                    	double bpy = ((double[])belief.get(1-i))[y1];
                        s += func.GetValue(y, y1) * bpy;
                    }
                    msg[y] = s;
                }
                NormalizeMessage();
            }
            ((Node) neighbor.get(i)).GetMessageFrom(nodeid, msg, gl);
        }
    }
    
    public void MaxSumPropagation(globalvar gl, boolean labeled_given)
    {
        for (int i = 0; i < 2; i ++)    
        {
            if (labeled_given &&((VariableNode) neighbor.get(i)).label_type == 0)
            {
                for (int y = 0; y < num_label; y ++)
                    msg[y] = 0;
                msg[((VariableNode) neighbor.get(i)).y] = 1.0;
            }
            else
            {
                for (int y = 0; y < num_label; y ++)
                {
                	
                	double bpy = ((double[])belief.get(1-i))[0];
                    double s = func.GetValue(y, 0) * bpy;
                    double tmp;
                    for (int y1 = 0; y1 < num_label; y1 ++)
                    {
                        bpy = ((double[])belief.get(1-i))[y1];
                        tmp = func.GetValue(y, y1) *bpy;
                        if (tmp > s) s = tmp;
                    }
                    msg[y] = s;
                }
                NormalizeMessage();
            }

            ((Node) neighbor.get(i)).GetMessageFrom(nodeid, msg, gl);
        }
    }
   
}
