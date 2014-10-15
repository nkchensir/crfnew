package crf;

import java.util.ArrayList;

public class DataNode
{

   public  int  label_type;
   public int   label;
   public int    num_attrib;
   public ArrayList<Integer> attrib;
   public ArrayList<Double>   value;

   public DataNode()
   {
	   label_type = 0;	
	   label =0;
	   num_attrib =0;
	   attrib = new ArrayList<Integer>();
	   value = new  ArrayList<Double>();
   }
   public int GetSize() { return 16 * (3 + num_attrib) + 16 * num_attrib; }
};