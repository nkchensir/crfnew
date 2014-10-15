package crf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

public class FactorGraph
{   
    public int    n, m, num_label;
    public int     num_node;

    public boolean               converged;
    public boolean                labeled_given;

    public VariableNode[]       var_node;
    public FactorNode[]         factor_node;
  //  public Node[][]              p_node;
    public ArrayList p_node;
    public ArrayList             bfs_node;

    // For each subgraph (connected component), we select one node as entry
    public ArrayList      entry;

    public int                 factor_node_used;

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
    
    public void InitGraph(int n, int m, int num_label)
    {
        labeled_given = false;
        this.n = n;
        this.m = m;
        this.num_label = num_label;
        this.num_node = n + m;
        
        var_node = new VariableNode[n];
        factor_node = new FactorNode[m];

        entry = new ArrayList();
        bfs_node = new ArrayList();
        int p_node_id = 0;

      //  p_node = new Node[n+m][];
        p_node = new ArrayList(); 
        for (int i = 0; i < n; i ++)
        {
        	var_node[i] = new VariableNode();
            var_node[i].nodeid = p_node_id;
            p_node_id++;
            p_node.add(var_node[i]);
            
            var_node[i].Init( num_label );
        }

        for (int i = 0; i < m; i ++)
        {
        	 factor_node[i] = new FactorNode();
            factor_node[i].nodeid = p_node_id;
            p_node_id++;
            p_node.add(factor_node[i]);
            factor_node[i].Init( num_label );
        }

        factor_node_used = 0;
    }
    public void AddEdge(int a, int b, EdgeFactorFunction func)
    {
        // AddEdge can be called at most m times
        if (factor_node_used == m) return;

        factor_node[factor_node_used].func = func;
            
        factor_node[factor_node_used].AddNeighbor( var_node[a] );
        factor_node[factor_node_used].AddNeighbor( var_node[b] );

        var_node[a].AddNeighbor( factor_node[factor_node_used] );
        var_node[b].AddNeighbor( factor_node[factor_node_used] );

        factor_node_used ++;
    }    
    public void GenPropagateOrder()
    {
        boolean[] mark = new boolean[num_node];
        bfs_node = new ArrayList() ;


        for (int i = 0; i < num_node; i ++)
            mark[i] = false;

        int head = 0, tail = -1;
        for (int i = 0; i < num_node; i ++)
        {
            if (! mark[i])
            {
                entry.add( p_node.get(i) );
                tail++;
                bfs_node.add(p_node.get(i));
                int type = ((Node)  p_node.get(i)).nodetype;
                mark[  ((Node)  p_node.get(i)).nodeid ] = true;

                while (head <= tail)
                {
                    Node u = (Node) bfs_node.get(head);
                    head++;
                    int ns = u.neighbor.size();
                    for (int s = 0;s < ns;s++)
                    if (! mark[((Node) u.neighbor.get(s)).nodeid])
                    {
                    	bfs_node.add(u.neighbor.get(s));
                    	tail++;
                    	mark[((Node) u.neighbor.get(s)).nodeid]=true;
                    }           
                    
                }
            }
        }
    }
    
    public void ClearDataForSumProduct()
    {   
        for (int i = 0; i < n; i ++)
        {
           DoubleArrFill(var_node[i].state_factor, num_label, 1.0 / num_label);          
        }

        for (int i = 0; i < num_node; i ++)
        {
            for (int t = 0; t < ((Node) p_node.get(i)).neighbor.size(); t ++)
            {
            	int type = ((Node) p_node.get(i)).nodetype;
            	if(type==1)
            		DoubleArrFill( (double[])((FactorNode) p_node.get(i)).belief.get(t), num_label, 1.0 / num_label);
            	else
            		DoubleArrFill( (double[])((VariableNode) p_node.get(i)).belief.get(t), num_label, 1.0 / num_label);	
            }
        }
    }
    public void ClearDataForMaxSum()
    {
        for (int i = 0; i < n; i ++)
        {
            DoubleArrFill(var_node[i].state_factor, num_label, 1.0 / num_label);
        }
        for (int i = 0; i < num_node; i ++)
        {
            for (int t = 0; t <((Node) p_node.get(i)).neighbor.size(); t ++)
            {
            	int type = ((Node) p_node.get(i)).nodetype;
            	if(type==1)
                for (int y = 0; y < num_label; y ++)
                	 ((double[])((FactorNode) p_node.get(i)).belief.get(t))[y]= 1.0 / num_label;
            	else
                    for (int y = 0; y < num_label; y ++)
            		 ((double[])((VariableNode) p_node.get(i)).belief.get(t))[y]= 1.0 / num_label;
            }
        }
    }
    
    public void SetVariableLabel(int u, int y) { var_node[u].y = y; }
    public void SetVariableStateFactor(int u, int y, double v){ var_node[u].state_factor[y] = v; }
    
    // Sum-Product
    public void BeliefPropagation(int max_iter,globalvar gl)
    {    
        int start, end, dir;

        converged = false;
        for (int iter = 0; iter < max_iter; iter ++)
        {
            gl.diff_max = 0.0;

            if (iter % 2 == 0)
            {
                start = num_node - 1;
                end = -1;
                dir = -1;
            }
            else
            {
                start = 0;
                end = num_node;
                dir = +1;
            }

            for (int p = start; p != end; p += dir)
            {
            	int type = ((Node) bfs_node.get(p)).nodetype;
            	if(type==1)
            		((FactorNode)bfs_node.get(p)).BeliefPropagation(gl, labeled_given);
            	else
            		 ((VariableNode)bfs_node.get(p)).BeliefPropagation(gl, labeled_given);	
            }

            if (gl.diff_max < 1e-6) break;
        }
    }
    public void CalculateMarginal()
    {
        for (int i = 0; i < n; i ++)
        {
            double sum_py = 0.0;
            for (int y = 0; y < num_label; y ++)
            {
              //  System.out.println(var_node[i].state_factor[y]);
                var_node[i].marginal[y] = var_node[i].state_factor[y];
                for (int t = 0; t < var_node[i].neighbor.size(); t ++)
                {
                    var_node[i].marginal[y] *= ((double[])var_node[i].belief.get(t))[y];
                   // System.out.println(((double[])var_node[i].belief.get(t))[y]);
                }
                sum_py += var_node[i].marginal[y];
            }
          //  System.out.println(sum_py);
            for (int y = 0; y < num_label; y ++)
            {
                var_node[i].marginal[y] /= sum_py;
            }
        }

        for (int i = 0; i < m; i ++)
        {
            double sump = 0.0;
            for (int a = 0; a < num_label; a ++)
                for (int b = 0; b < num_label; b ++)
                {
                    factor_node[i].marginal[a][b] +=
                    		((double[])factor_node[i].belief.get(0))[a]
                          * ((double[])factor_node[i].belief.get(1))[b]
                          * factor_node[i].func.GetValue(a, b);
                    sump += factor_node[i].marginal[a][b];
                }
            for (int a = 0; a < num_label; a ++)
                for (int b = 0; b < num_label; b ++)
                    factor_node[i].marginal[a][b] /= sump;
        }
    }

    // Max-Sum
    public void MaxSumPropagation(int max_iter,globalvar gl)
    {
        int start, end, dir;

        converged = false;
        for (int iter = 0; iter < max_iter; iter ++)
        {
            gl.diff_max = 0;

            if (iter % 2 == 0)
            {
                start = num_node - 1;
                end = -1;
            	dir = -1;
            }
            else
            {
                start = 0;
                end = num_node;
                dir = +1;
            }

            for (int p = start; p != end; p += dir)
            {                
            	int type = ((Node) bfs_node.get(p)).nodetype;
            	if(type==1)
            		((FactorNode)bfs_node.get(p)).MaxSumPropagation(gl, labeled_given);
            	else
            		 ((VariableNode)bfs_node.get(p)).MaxSumPropagation(gl, labeled_given);	
            }

            if (gl.diff_max < 1e-6) break;
        }
    }
   
};