package crf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import crf.Util.*;

public class Node {

	public int nodeid;
	public int num_label;
    public int nodetype;
    public ArrayList   neighbor;
    public ArrayList   belief;
    public  HashMap<Integer,Integer>  neighbor_pos;
	public double[] msg;

	public Node()
	{
	//msg = new double[];
		 nodeid=-1;
		 num_label=-1;
		 nodetype=-1;
		 neighbor = new ArrayList();
		 belief = new ArrayList();
		 neighbor_pos = new  HashMap<Integer,Integer>();
		 
	}
	public void  NodeInit()
	{
	//msg = new double[];
		 nodeid=-1;
		 num_label=-1;
		 nodetype=-1;
		 neighbor = new ArrayList();
		 belief = new ArrayList();
		 neighbor_pos = new  HashMap<Integer,Integer>();
		 
	}
    public double[] GetDoubleArr(int size)
    {
        double[] arr = new double[size];
        return arr;
    }
    
    public void DoubleArrFill(double[] arr, int size, double v)
    {
        for (int i = 0; i < size; i ++)
            arr[i] = v;
    }
    
	public void  BasicInit(int num_label)
	{
	    this.num_label = num_label;
	    msg = new double[num_label];
	    neighbor = new ArrayList();
	    belief = new ArrayList();
	}

	public void NormalizeMessage()
	{
	    double s = 0.0;
	    for (int y = 0; y < num_label; y ++)
	        s += msg[y];
	    for (int y = 0; y < num_label; y ++)
	        msg[y] /= s;
	}

	
	public void AddNeighbor(Node ng)
	 {
	        neighbor_pos.put(ng.nodeid, neighbor.size());
	        neighbor.add(ng);

	        belief.add( GetDoubleArr(num_label) );
	 }


	public void BeliefPropagation(double diff_max, boolean labeled_given){}
	public void MaxSumPropagation(double diff_max, boolean labeled_given){}

	    void GetMessageFrom(int u, double[] msgvec, globalvar gl)
	    {
	        int p = neighbor_pos.get(u);
	        for (int y = 0; y < num_label; y ++)
	        {
	        	double bpy = ((double[])belief.get(p))[y];
	            if (Math.abs( bpy- msgvec[y]) > gl.diff_max)
	                 gl.diff_max = Math.abs(bpy - msgvec[y]);
	            ((double[])belief.get(p))[y] = msgvec[y];
	        }
	    }
	    
}
