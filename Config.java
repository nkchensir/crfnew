package crf;

public class Config
{
    public int   my_rank;
    public int   num_procs;

    public String      task;

    public String      train_file;
    public String      test_file;
    public String      pred_file;

    public String      dict_file;
    public String      src_model_file;
    public String      dst_model_file;

    public double      eps;

    public int         max_iter;
    public int         max_bp_iter;

    public double      gradient_step;

    public boolean        has_attrib_value;
    public int         optimization_method;

    public boolean        eval_each_iter;

    public double      penalty_sigma_square;
    

    Config(){ SetDefault(); }
    public void SetDefault()
    {
        max_iter = 1000;
        max_bp_iter = 20;
        this.my_rank =0;
        this.task = "-est";
        this.train_file = "train.txt";
        this.test_file = "test.txt";

        this.dict_file = "dict.txt";
        this.pred_file = "pred.txt";

        //this.train_file = "scene_data\\train.txt";
        //this.test_file = "scene_data\\test.txt";

        this.train_file = "/home/cwz/GCRF/src/data/name3crf.txt";
//this.test_file = "zz_data\\run\\test.txt";
        
        //this.train_file = "relation_graph\\1ha";
        //this.test_file = "relation_graph\\test.txt";

        //this.train_file = "relation_ind\\train.txt";
        //this.test_file = "relation_ind\\test.txt";

        //this.train_file = "new_adv_data\\parta.dat";
        //this.test_file  = "new_adv_data\\partb.dat";
        
        //this.has_attrib_value = false;
        this.has_attrib_value = true;
        this.eps = 1e-3;

        //this.optimization_method = LBFGS;
        this.optimization_method = 1;
        this.gradient_step =1;

        this.dst_model_file = "model-final.txt";

        this.eval_each_iter = true;

        this.penalty_sigma_square = 0.00001;
    }

    // false => parameter wrong
  //  public boolean LoadConfig(int my_rank, int num_procs, int argc, char* argv[]);
   // public void ShowUsage();
};
