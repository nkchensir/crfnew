package crf;

import java.io.IOException;

public class OpenCRFMain {

	public void Estimate(Config conf) throws IOException
	{
	    GlobalDataSet dataset;
	   // DataSet* dataset;

	       dataset = new GlobalDataSet();
	       dataset.LoadData(conf.train_file, conf);
	       dataset.label_dict.SaveMappingDict(conf.dict_file);

	   System.out.println("num_label = "+Integer.toString(dataset.num_label));
	  // System.out.println("num_sample = "+Integer.toString(dataset.num_sample));
	   System.out.println("num_edge_type = "+Integer.toString(dataset.num_edge_type));
	   System.out.println("num_attrib_type ="+Integer.toString(dataset.num_attrib_type));
	   
	   dataset.num_node = dataset.node.size();
	   dataset.num_edge = dataset.edge.size();
	   System.out.println("num node ="+Integer.toString(dataset.num_node));
	   System.out.println("num edge ="+Integer.toString(dataset.num_edge));
	    CRFModel model = new CRFModel();
	    
	    model.InitTrain(conf, dataset);    

	    System.out.println("Start Training...\n");
		//flush(stdout);
	    model.Train();
	    
	  //  if (conf.my_rank == 0)
	  ////    model.SaveModel(conf.dst_model_file);
	      
	    //MakeEvaluate(conf, g_dataset, model);
	}
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
	    Config conf = new Config();
	    OpenCRFMain  My = new OpenCRFMain();
        My.Estimate(conf);
	}

}
