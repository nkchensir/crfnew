package crf;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class CRFModel
{

    public Config     conf;
	public GlobalDataSet   train_data;
	public GlobalDataSet    test_data;

    public int         num_sample;
    public int         num_label;
    public int         num_attrib_type;
    public int         num_edge_type;
            
    public int         num_feature;
    public globalvar gl;
    public double[]      lambda;
    public FactorGraph sample_factor_graph;

    public int             num_attrib_parameter;
    public int             num_edge_feature_each_type;
    public HashMap<Integer, Integer>   edge_feature_offset;
    public EdgeFactorFunction[][]   func_list;
    public EdgeFactorFunction   func;

//    public CRFModel(){}

    public void InitTrain(Config  conf, GlobalDataSet train_data)
    {
        this.conf = new Config();
       // this.conf = conf;
        this.train_data = train_data;

        /*
        public int   num_node;
        public int    num_edge;
        public ArrayList<DataNode>   node;
        public ArrayList<DataEdge> edge;
        
        int num_label;
        int num_attrib_type;
        int num_edge_type;*/
        gl = new globalvar();
        edge_feature_offset = new HashMap<Integer, Integer>();
        num_sample = 1;
        num_label = train_data.num_label;
        num_attrib_type = train_data.num_attrib_type;
        num_edge_type = train_data.num_edge_type;
        num_edge_type = 1;
        GenFeature();
        lambda = new double[num_feature];
        // Initialize parameters
        for (int i = 0; i < num_feature; i ++)
            lambda[i] = 0.0;
        SetupFactorGraphs();
    }
    public void GenFeature()
    {
        num_feature = 0;

        // state feature: f(y, x)
        num_attrib_parameter = num_label * num_attrib_type;
        num_feature += num_attrib_parameter;

        // edge feature: f(edge_type, y1, y2)
        edge_feature_offset.clear();
        int offset = 0;
        for (int y1 = 0; y1 < num_label; y1 ++)
            for (int y2 = y1; y2 < num_label; y2 ++)
            {
                edge_feature_offset.put(y1 * num_label + y2, offset);
                offset ++;
            }
        num_edge_feature_each_type = offset;
        num_feature += num_edge_type * num_edge_feature_each_type;
    }
    public void SetupFactorGraphs()
    {
        double[]  p_lambda =Arrays.copyOfRange(lambda,num_attrib_parameter,num_attrib_parameter+ num_edge_feature_each_type+1);
      //  func_list = new EdgeFactorFunction*[ num_edge_type ];
        func = new EdgeFactorFunction(num_label, p_lambda, edge_feature_offset);
        
        /*
        for (int i = 0; i < num_edge_type; i ++)
        {
            func_list[i] = new EdgeFactorFunction(num_label, p_lambda, &edge_feature_offset);
            p_lambda += num_edge_feature_each_type;
        }*/

        sample_factor_graph = new FactorGraph();
        	int s=0;

            int n = train_data.num_node;
            int m = train_data.num_edge;
            
            sample_factor_graph.InitGraph(n, m, num_label);

            // Add node info
            for (int i = 0; i < n; i ++)
            {
                sample_factor_graph.SetVariableLabel(i, (train_data.node).get(i).label);
                sample_factor_graph.var_node[i].label_type = (train_data.node).get(i).label_type;
            }   
            
            // Add edge info
            for (int i = 0; i < m; i ++)
            {
                sample_factor_graph.AddEdge((train_data.edge).get(i).a, (train_data.edge).get(i).b, func);
            }

            sample_factor_graph.GenPropagateOrder();
        
    }

    public void Train() throws IOException
    {    
        double[] gradient;
        double  f;          // log-likelihood

        gradient = new double[num_feature + 1];

        // Master node
    //    if (conf.my_rank == 0)
        {
            ///// Initilize all info

            // Data Varible         
            double  old_f = 0.0;

            // Variable for lbfgs optimization
            int     m_correlation = 3;
            double[] work_space = new double[num_feature * (2 * m_correlation + 1) + 2 * m_correlation];
            int     diagco = 0;
            double[] diag = new double[num_feature];
            int[]     iprint = {-1, 0}; // do not print anything
            double  eps = 0;
            eps = conf.eps;
            double  xtol = 1.0e-16;
            int     iflag = 0;
            long startMili = 0;
            // Other Variables
            int     num_iter;
            double[] tmp_store = new double[num_feature + 1];
            
            // Main-loop of CRF
            // Paramater estimation via L-BFGS
            num_iter = 0;

            double start_time, end_time;

            do {
                num_iter ++;

              //  if (conf.my_rank == 0)
                {
                  startMili=System.currentTimeMillis();
                }
                            
                // Step A. Send lambda to all procs
             //   Transmitter::Master_SendDoubleArray(lambda, num_feature, conf.num_procs);

                // Step B. Calc gradient and log-likehood of the local datas
                f = CalcGradient(gradient);

                // Step C. Collect gradient and log-likehood from all procs
              //  Transmitter::Master_CollectGradientInfo(gradient, &f, num_feature, tmp_store, conf.num_procs);

                // Step 4. Opitmization by L-BFGS
               System.out.printf("[Iter %3d] log-likelihood : %f\n", num_iter, f);
            //    fflush(stdout);

                // If diff of log-likelihood is small enough, break.
                if (Math.abs(old_f - f) < eps) break;
                	old_f = f;
    		
    		    // Negate f and gradient vector because the LBFGS optimization below minimizes the ojective function while we would like to maximize it
                f *= -1;
                for (int i = 0; i < num_feature; i ++)
                    gradient[i] *= -1;

                // Invoke L-BFGS

                if (conf.optimization_method == 0)
                {
                    //lbfgs_(&num_feature, &m_correlation, lambda, &f, gradient, &diagco, diag, iprint, &eps, &xtol, work_space, &iflag);

                    // Checking after calling LBFGS
                    if (iflag < 0) // LBFGS error
    		        {
    		            System.out.printf("LBFGS routine encounters an error\n");
    		            break;
    		        }
                }
                else
                {
                    // Normalize Graident
                    double g_norm = 0.0;
                    for (int i = 0; i < num_feature; i ++)
                        g_norm += gradient[i] * gradient[i];
                    g_norm = Math.sqrt(g_norm);
                    
                    if (g_norm > 1e-8)
                    {
                        for (int i = 0; i < num_feature; i ++)
                            gradient[i] /= g_norm;
                    }

                    for (int i = 0; i < num_feature; i ++)
                        lambda[i] -= gradient[i] * conf.gradient_step;
                    iflag = 1;
                }

                if (conf.eval_each_iter && conf.my_rank == 0)
                {
                    //SelfEvaluate();
                }

           //     if (conf.my_rank == 0)
                {
                //    end_time = clock() / (double)CLK_TCK;
                    long endMili=System.currentTimeMillis();
                //    FILE* ftime = fopen("time.out", "a");
                //    fprintf(ftime, "start_time = %.6lf\n", start_time);
                //    fprintf(ftime, "end_time = %.6lf\n", end_time);
                //    fprintf(ftime, "cost = %.6lf\n", end_time - start_time);
                    
                //    fclose(ftime);

                   System.out.printf("!!! Time cost = %d\n", endMili - startMili);
                    //fflush(stdout);
                }
            } while (iflag != 0 && num_iter < conf.max_iter);



         //   Transmitter::Master_SendQuit(conf.num_procs);

         //   delete[] tmp_store;

         //   delete[] work_space;
        //    delete[] diag;
        }

     //   delete[] gradient;
    }

    public double CalcGradient(double[] gradient) throws IOException
    {
        double  f;
            
        // Initialize

        // If there is a square penalty, gradient should be initialized with (- lambda[i] / sigma^2). 
        // f should be accordingly modified as : -||lambda||^2/ (2*sigma^2)
        // note : should be added only in one procs (master)

        f = 0.0;
        for (int i = 0; i < num_feature; i ++)
        {
        //    gradient[i] = - lambda[i] / conf.penalty_sigma_square;
            gradient[i] = 0; // no penalty
        }

        // Calculation
       // System.out.println(gradient[0]);
        for (int i = 0; i < num_sample; i ++)
        {
            // double t = CalcGradientForSample(train_data.sample[i], &sample_factor_graph[i], gradient);
            double t = CalcPartialLabeledGradientForSample(train_data, sample_factor_graph, gradient);
            f += t;		 
        }
      //  System.out.println(gradient[0]);
        return f;
    }
    public double CalcGradientForSample(GlobalDataSet sample, FactorGraph factor_graph, double[] gradient)
{
    factor_graph.ClearDataForSumProduct();

    // Set state_factor
    int n = sample.num_node;
    int m = sample.num_edge;

    for (int i = 0; i < n; i ++)
    {
     //   double* p_lambda = lambda;
    	int pos = 0;
        for (int y = 0; y < num_label; y ++)
        {
            double v = 1;
            for (int t = 0; t < sample.node.get(i).num_attrib; t ++)
                v *= Math.exp( lambda[ sample.node.get(i).attrib.get(t) ] * sample.node.get(i).value.get(t) );
            factor_graph.SetVariableStateFactor(i, y, v);
            pos+= num_attrib_type;
          //  p_lambda += num_attrib_type;
        }
    }    

    factor_graph.BeliefPropagation(conf.max_bp_iter,gl);
    factor_graph.CalculateMarginal();
    
    // Calculate gradient & log-likelihood
    double f = 0.0, Z = 0.0;

    // \sum \lambda_i * f_i
    for (int i = 0; i < n; i ++)
    {
        int y = sample.node.get(i).label;
        for (int t = 0; t < sample.node.get(i).num_attrib; t ++)
            f += lambda[this.GetAttribParameterId(y, sample.node.get(i).attrib.get(t))] * sample.node.get(i).value.get(t);
    }
    for (int i = 0; i < m; i ++)
    {
        int a = sample.node.get(sample.edge.get(i).a).label;
        int b = sample.node.get(sample.edge.get(i).b).label;        
        f += lambda[this.GetEdgeParameterId(sample.edge.get(i).edge_type, a, b)];
    }

    // calc log-likelihood
    //  using Bethe Approximation
    for (int i = 0; i < n; i ++)
    {
        for (int y = 0; y < num_label; y ++)
        {
            for (int t = 0; t < sample.node.get(i).num_attrib; t ++)
                Z += lambda[this.GetAttribParameterId(y, sample.node.get(i).attrib.get(t))] * sample.node.get(i).value.get(t) * factor_graph.var_node[i].marginal[y];
        }
    }
    for (int i = 0; i < m; i ++)
    {
        for (int a = 0; a < num_label; a ++)
            for (int b = 0; b < num_label; b ++)
            {
                Z += lambda[this.GetEdgeParameterId(sample.edge.get(i).edge_type, a, b)] * factor_graph.factor_node[i].marginal[a][b];
            }
    }
    // Edge entropy
    for (int i = 0; i < m; i ++)
    {
        double h_e = 0.0;
        for (int a = 0; a < num_label; a ++)
            for (int b = 0; b < num_label; b ++)
            {
                if (factor_graph.factor_node[i].marginal[a][b] > 1e-10)
                    h_e += - factor_graph.factor_node[i].marginal[a][b] * Math.log(factor_graph.factor_node[i].marginal[a][b]);
            }
        Z += h_e;
    }
    // Node entroy
    for (int i = 0; i < n; i ++)
    {
        double h_v = 0.0;
        for (int a = 0; a < num_label; a ++)
            if (Math.abs(factor_graph.var_node[i].marginal[a]) > 1e-10)
                h_v += - factor_graph.var_node[i].marginal[a] * Math.log(factor_graph.var_node[i].marginal[a]);
        Z -= h_v * ((int)factor_graph.var_node[i].neighbor.size() - 1);
    }

    f -= Z;
 //   fflush(stdout);

    // calc gradient
    for (int i = 0; i < n; i ++)
        for (int t = 0; t < sample.node.get(i).num_attrib; t ++)
            gradient[ GetAttribParameterId(sample.node.get(i).label, sample.node.get(i).attrib.get(t)) ] += sample.node.get(i).value.get(t);
    for (int i = 0; i < m; i ++)
        gradient[ GetEdgeParameterId(sample.edge.get(i).edge_type, 
                                     sample.node.get(sample.edge.get(i).a).label, 
                                     sample.node.get(sample.edge.get(i).b).label)] += 1.0;

    // - expectation
    for (int i = 0; i < n; i ++)
    {
        for (int y = 0; y < num_label; y ++)
        {
            for (int t = 0; t < sample.node.get(i).num_attrib; t ++)
                gradient[ GetAttribParameterId(y, sample.node.get(i).attrib.get(t)) ] -= sample.node.get(i).value.get(t) * factor_graph.var_node[i].marginal[y];
        }
    }
    for (int i = 0; i < m; i ++)
    {
        for (int a = 0; a < num_label; a ++)
            for (int b = 0; b < num_label; b ++)
            {
                gradient[ GetEdgeParameterId(sample.edge.get(i).edge_type, a, b) ] -= factor_graph.factor_node[i].marginal[a][b];
            }
    }

    return f;
}
   
    public double CalcPartialLabeledGradientForSample(GlobalDataSet sample, FactorGraph factor_graph, double[] gradient) throws IOException
    {   
        int n = sample.num_node;
        int m = sample.num_edge;

        //****************************************************************
        // Belief Propagation 1: labeled data are given.
        //****************************************************************

        factor_graph.labeled_given = true;
        factor_graph.ClearDataForSumProduct();

        // Set state_factor
        for (int i = 0; i < n; i ++)
        {
          //  double[] p_lambda = lambda;
            int pos = 0;
            for (int y = 0; y < num_label; y ++)
            {
                if ((sample.node).get(i).label_type == 0 && y != (sample.node).get(i).label)
                {
                    factor_graph.SetVariableStateFactor(i, y, 0);
                }
                else
                {
                    double v = 1;
                    for (int t = 0; t < (sample.node).get(i).num_attrib; t ++)
                        v *= Math.exp( lambda[ pos+((sample.node).get(i).attrib).get(t) ] * ((sample.node).get(i).value).get(t) );   
                    	//System.out.println("v="+Double.toString(v));
                    factor_graph.SetVariableStateFactor(i, y, v);
                }
                pos += num_attrib_type;
                //p_lambda += num_attrib_type;
              //  p_lambda = Arrays.copyOfRange(lambda,num_attrib_type)
            }
        }

        factor_graph.BeliefPropagation(conf.max_bp_iter,gl);
        factor_graph.CalculateMarginal();    

        /***
         * Gradient = E_{Y|Y_L} f_i - E_{Y} f_i
         */

        // calc gradient part : + E_{Y|Y_L} f_i
        for (int i = 0; i < n; i ++)
        {
            for (int y = 0; y < num_label; y ++)
            {
                for (int t = 0; t < (sample.node).get(i).num_attrib; t ++)
                {
                    gradient[ GetAttribParameterId(y, sample.node.get(i).attrib.get(t)) ] += sample.node.get(i).value.get(t) * factor_graph.var_node[i].marginal[y];
               //     System.out.println("v="+Double.toString(gradient[ GetAttribParameterId(y, sample.node.get(i).attrib.get(t)) ]));
                }
            }
        }

        for (int i = 0; i < m; i ++)
        {
            for (int a = 0; a < num_label; a ++)
                for (int b = 0; b < num_label; b ++)
                {
                    gradient[ GetEdgeParameterId(sample.edge.get(i).edge_type, a, b) ] += factor_graph.factor_node[i].marginal[a][b];
                  
                }
        }

        //****************************************************************
        // Belief Propagation 2: labeled data are not given.
        //****************************************************************


        factor_graph.ClearDataForSumProduct();
        factor_graph.labeled_given = false;

        for (int i = 0; i < n; i ++)
        {
           // double* p_lambda = lambda;
        	int pos = 0;
            for (int y = 0; y < num_label; y ++)
            {
                double v = 1;
                for (int t = 0; t < sample.node.get(i).num_attrib; t ++)
                    v *= Math.exp( lambda[pos+ sample.node.get(i).attrib.get(t) ] * sample.node.get(i).value.get(t) );
             //   System.out.println("v="+Double.toString(v));
                factor_graph.SetVariableStateFactor(i, y, v);
                pos += num_attrib_type;
            }
        }    

        factor_graph.BeliefPropagation(conf.max_bp_iter,gl);
        factor_graph.CalculateMarginal();
            
        // calc gradient part : - E_{Y} f_i
        for (int i = 0; i < n; i ++)
        {
            for (int y = 0; y < num_label; y ++)
            {
                for (int t = 0; t < sample.node.get(i).num_attrib; t ++)
                    gradient[ GetAttribParameterId(y, sample.node.get(i).attrib.get(t)) ] -= sample.node.get(i).value.get(t) * factor_graph.var_node[i].marginal[y];
            }
        }
        for (int i = 0; i < m; i ++)
        {
            for (int a = 0; a < num_label; a ++)
                for (int b = 0; b < num_label; b ++)
                {
                    gradient[ GetEdgeParameterId(sample.edge.get(i).edge_type, a, b) ] -= factor_graph.factor_node[i].marginal[a][b];
                }
        }
        
        // Calculate gradient & log-likelihood
        double f = 0.0, Z = 0.0;

        // \sum \lambda_i * f_i
        for (int i = 0; i < n; i ++)
        {
            int y = sample.node.get(i).label;
            for (int t = 0; t < sample.node.get(i).num_attrib; t ++)
                f += lambda[this.GetAttribParameterId(y, sample.node.get(i).attrib.get(t))] * sample.node.get(i).value.get(t);
        }
        for (int i = 0; i < m; i ++)
        {
            int a = sample.node.get( sample.edge.get(i).a ).label;
            int b = sample.node.get( sample.edge.get(i).b ).label;
            f += lambda[this.GetEdgeParameterId(sample.edge.get(i).edge_type, a, b)];
        }
      //  System.out.println("f="+Double.toString(f));
        // calc log-likelihood
        //  using Bethe Approximation
        for (int i = 0; i < n; i ++)
        {
            for (int y = 0; y < num_label; y ++)
            {
                for (int t = 0; t < sample.node.get(i).num_attrib; t ++)
                    Z += lambda[this.GetAttribParameterId(y, sample.node.get(i).attrib.get(t))] * sample.node.get(i).value.get(t) * factor_graph.var_node[i].marginal[y];
            }
        }
      //  System.out.println("Z="+Double.toString(Z));
        for (int i = 0; i < m; i ++)
        {
            for (int a = 0; a < num_label; a ++)
                for (int b = 0; b < num_label; b ++)
                {
                    Z += lambda[this.GetEdgeParameterId(sample.edge.get(i).edge_type, a, b)] * factor_graph.factor_node[i].marginal[a][b];
                }
        }
     //   System.out.println("Z="+Double.toString(Z));
        // Edge entropy
        for (int i = 0; i < m; i ++)
        {
            double h_e = 0.0;
            for (int a = 0; a < num_label; a ++)
                for (int b = 0; b < num_label; b ++)
                {
                    if (factor_graph.factor_node[i].marginal[a][b] > 1e-10)
                        h_e += - factor_graph.factor_node[i].marginal[a][b] * Math.log(factor_graph.factor_node[i].marginal[a][b]);
                //    System.out.printf("fm=%f\n",factor_graph.factor_node[i].marginal[a][b]);
                }
            Z += h_e;
        }
     //   System.out.println("log2.7="+Double.toString( Math.log(2.7)));
      //  System.out.println("Z="+Double.toString(Z));
        // Node entroy
        for (int i = 0; i < n; i ++)
        {
            double h_v = 0.0;
            for (int a = 0; a < num_label; a ++)
                if (Math.abs(factor_graph.var_node[i].marginal[a]) > 1e-10)
                    h_v += - factor_graph.var_node[i].marginal[a] * Math.log(factor_graph.var_node[i].marginal[a]);
            Z -= h_v * ((int)factor_graph.var_node[i].neighbor.size() - 1);
        }
     //   System.out.println("Z="+Double.toString(Z));
        f -= Z;
        
    //#ifdef DO_EVAL
        // Let's take a look of current accuracy

        factor_graph.ClearDataForMaxSum();
        factor_graph.labeled_given = true;

            for (int i = 0; i < n; i ++)
            {
             //   double* p_lambda = lambda;
            	int pos = 0;
                for (int y = 0; y < num_label; y ++)
                {
                    double v = 1.0;
                    for (int t = 0; t < sample.node.get(i).num_attrib; t ++)
                        v *=Math.exp( lambda[ pos+sample.node.get(i).attrib.get(t) ] * sample.node.get(i).value.get(t) );
                    factor_graph.SetVariableStateFactor(i, y, v);
                    
                    pos += num_attrib_type;
                 //   p_lambda += num_attrib_type;
                }
            }    

            factor_graph.MaxSumPropagation(conf.max_bp_iter,gl);

            int[] inf_label = new int[n];
    		double[][] label_prob = new double[num_label][];
    		for (int p = 0; p < num_label; p ++)
    			label_prob[p] = new double[n];

            for (int i = 0; i < n; i ++)
            {
                int ybest = -1;
                double vbest = -1, v;
    			double vsum = 0.0;
                for (int y = 0; y < num_label; y ++)
                {
                    v = factor_graph.var_node[i].state_factor[y];
                    for (int t = 0; t < factor_graph.var_node[i].neighbor.size(); t ++)
                        v *= ((double[]) factor_graph.var_node[i].belief.get(t))[y];
                    if (ybest < 0 || v > vbest)
                    {
                        ybest = y;
                        vbest = v;
                    }

    				label_prob[y][i] = v;
    				vsum += v;
                }

                inf_label[i] = ybest;

    			for (int y = 0; y < num_label; y ++)
    				label_prob[y][i] /= vsum;
            }

        int hit = 0, miss = 0;
        int hitu = 0, missu = 0;

    	int[][] cnt = new int[10][10];
    	int[][] ucnt = new int[10][10];


    	
		BufferedWriter out = new BufferedWriter(new FileWriter( new File("pred.txt")));
		
		/*
        for (int i = 0; i < keys.size(); i ++)
            out.write(keys.get(i)+Integer.toString(i)+"\n");			
		out.flush();
		out.close();
		*/
    //	FILE *pred_out = fopen("pred.txt", "w");

        for (int i = 0; i < n; i ++)
        {
            /*
            int ymax = 0;
            for (int y = 0; y < num_label; y ++)
                if (factor_graph.var_node[i].marginal[y] > factor_graph.var_node[i].marginal[ymax])
                    ymax = y;
                    */
        	  out.write(Integer.toString(inf_label[i])+" "+Double.toString(label_prob[inf_label[i]][i])+"\n");			
    	//	fprintf(pred_out, "%d %f\n", inf_label[i],label_prob[inf_label[i]][i]);
            if (inf_label[i] == sample.node.get(i).label)
                hit ++;
            else
                miss ++;

    		cnt[ inf_label[i] ][ sample.node.get(i).label ] ++;

            if (sample.node.get(i).label_type == 1)
            {
                if (inf_label[i] == sample.node.get(i).label)
                    hitu ++;
                else
                    missu ++;

    			ucnt[ inf_label[i] ][ sample.node.get(i).label ] ++;
            }
        }
		out.close();
        //printf("HIT = %4d, MISS = %4d, All_Accuracy = %.5lf Unknown_Accuracy = %.5lf\n", hit, miss, (double)hit / (hit + miss), (double)hitu / (hitu + missu));

 
        {
        	System.out.printf("train precision = %4.5f\n", (double)(hit -hitu)/(miss+ hit-hitu-missu));
        	System.out.printf("test precision = %4.5f\n", (double)(hitu)/(hitu+missu));
        	/*
            int[] dat = new int[12];
            //memset(dat, 0, sizeof(dat));
            //Transmitter::Master_CollectIntArr(dat, 12, conf.num_procs);

            hit += dat[0]; hitu += dat[1];
            miss += dat[2]; missu += dat[3];
            cnt[0][0] += dat[4]; cnt[0][1] += dat[5]; cnt[1][0] += dat[6]; cnt[1][1] += dat[7];
            ucnt[0][0] += dat[8]; ucnt[0][1] += dat[9]; ucnt[1][0] += dat[10]; ucnt[1][1] += dat[11];

    	    System.out.println("A_HIT  = "+Integer.toString(hit)+"  U_HIT  = "+Integer.toString(hitu));
    	    System.out.println("A_MISS = "+Integer.toString(miss)+" miss U_MISS = \n"+ Integer.toString(missu));

    	    //!!!!!!!! make sure, the first instance is "positive"

    	    // 0 . positive
    	    // 1 . negative

    	    double ap = (double)cnt[0][0] / (cnt[0][0] + cnt[0][1]);
    	    double up = (double)ucnt[0][0] / (ucnt[0][0] + ucnt[0][1]);

    	    double ar = (double)cnt[0][0] / (cnt[0][0] + cnt[1][0]);
    	    double ur = (double)ucnt[0][0] / (ucnt[0][0] + ucnt[1][0]);

    	    double af = 2 * ap * ar / (ap + ar);
    	    double uf = 2 * up * ur / (up + ur);

    	   System.out.printf("A_Accuracy  = %4.5f     U_Accuracy  = %4.5f\n", (double)hit / (hit + miss), (double)hitu / (hitu + missu));
    	   System.out.printf("A_Precision = %4.5f     U_Precision = %4.5f\n", ap, up);
    	   System.out.printf("A_Recall    = %4.5f     U_Recall    = %4.5f\n", ar, ur);
    	   System.out.printf("A_F1        = %4.5f     U_F1        = %4.5f\n", af, uf);
    		*/
          //  fflush(stdout);

        }

         out = new BufferedWriter(new FileWriter( new File("uncertainty.txt")));
        
    //	FILE* fprob = fopen("uncertainty.txt", "w");
    	for (int i = 0; i < n; i ++)
    	{
    		if (sample.node.get(i).label_type == 0)
    		{
    			for (int y = 0; y < num_label; y ++)
    				out.write("-1 ");
    			out.write("\n");
    		}
    		else
    		{
    			for (int y = 0; y < num_label; y ++)
    				out.write(Double.toString(label_prob[y][i]));
    			//	fprintf(fprob, "%4.2f ", label_prob[y][i]);
    			out.write("\n");
    		}
    	}
    	out.close();
    //	fclose(fprob);
        

    	//delete[] inf_label;
    //	for (int y = 0; y < num_label; y ++)
    //		delete[] label_prob[y];
    //	delete[] label_prob;
    //#endif

        return f;
    }    

    //void SelfEvaluate();
   // void PartialSelfEvaluation();
    
  //  void InitEvaluate(Config* conf, DataSet* test_data);
  //  void Evalute();

    int GetAttribParameterId(int y, int x){ return y * num_attrib_type + x; }
    int GetEdgeParameterId(int edge_type, int a, int b)
    { 
        int offset = edge_feature_offset.get((a<b?a:b) * num_label + (a>b?a:b));
        return num_attrib_parameter + edge_type * num_edge_feature_each_type + offset;
    }

   // ~CRFModel() { Clean(); }
  //  void Clean();

  //  void SaveModel(const char* file_name);
 //   void LoadModel(const char* file_name);
};