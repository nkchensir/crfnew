package crf;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class GlobalDataSet
{
 //   vector<DataSample*> sample;

    public int   num_node;
    public int    num_edge;
    public ArrayList<DataNode>   node;
    public ArrayList<DataEdge> edge;
    
    int num_label;
    int num_attrib_type;
    int num_edge_type;

    MappingDict         label_dict;
    MappingDict         attrib_dict;
    MappingDict         edge_type_dict;

    public GlobalDataSet()
    {
    	node = new ArrayList<DataNode>();
    	edge = new ArrayList<DataEdge>();
    	label_dict = new MappingDict();
    	attrib_dict = new MappingDict();
    	 edge_type_dict = new MappingDict();
    }
   public  void LoadData(String data_file, Config conf) throws IOException
   {
	   

       FileReader fr=new FileReader(data_file);
       BufferedReader br=new BufferedReader(fr);
       String line="";
       String[] tokens=null;
       while ((line=br.readLine())!=null) 
       {
    	  // System.out.println("newline: "+line.substring(0,100));
           tokens=line.split(" ");

	        if (tokens[0].equals("#edge")) //edge
	        {
	          //  System.out.println("tokens[0]:  "+tokens[0]);
	            DataEdge  curt_edge = new DataEdge();

	            curt_edge.a = Integer.parseInt(tokens[1]) - 1;
	            curt_edge.b = Integer.parseInt(tokens[2]) - 1;

	            edge.add(curt_edge);
	            //curt_edge->edge_type = 0;
	        //    if (tokens.size() >= 4)
	          //      curt_edge->edge_type = edge_type_dict.GetId( tokens[3] );//edgetype都不要了

	         //   curt_sample->edge.push_back(curt_edge);
	        }
	        else
	        {
	            DataNode curt_node = new DataNode();

	            char label_type = tokens[0].charAt(0);
	            String label_name = tokens[0].substring(1);

	            curt_node.label = label_dict.GetId(label_name);
	            if (label_type == '+')
	                curt_node.label_type = 0;
	            else if (label_type == '?')
	                curt_node.label_type = 1;
	            else {
	                System.out.println("Data format wrong! Label must start with +/?\n");

	            }
	            
	            for (int i = 1; i < tokens.length; i ++)
	            {
	                if (conf.has_attrib_value)//default true
	                {

	                    String[] key_val =null;	
	                    key_val = tokens[i].split(":");
	                    curt_node.attrib.add( attrib_dict.GetId(key_val[0]) );
	                    curt_node.value.add( Double.parseDouble(key_val[1]) );
	                }
	                else
	                {
	                    curt_node.attrib.add( attrib_dict.GetId(tokens[i]) );
	                    curt_node.value.add(1.0);
	                }
	            }

	            curt_node.num_attrib = curt_node.attrib.size();
	            node.add( curt_node );
	        }
           
        //   System.out.println(arrs[0] + " : " + arrs[1] + " : " + arrs[2]);
       }
       
	    num_label = label_dict.GetSize();
	    num_attrib_type = attrib_dict.GetSize();
	    num_edge_type = edge_type_dict.GetSize();
	    if (num_edge_type == 0) num_edge_type = 1;

       br.close();
       fr.close();
       

	}

 //   void LoadDataWithDict(const char* data_file, Config* conf, const MappingDict& ref_label_dict, const MappingDict& ref_attrib_dict, const MappingDict& ref_edge_type_dict);
};
