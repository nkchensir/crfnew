package crf;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Vector;
import java.io.IOException;
import javax.swing.JComboBox.KeySelectionManager;

public class MappingDict
{
    HashMap<String, Integer>    dict;
    Vector<String>      keys;
    
    public MappingDict()
    {
    dict = new HashMap<String, Integer>();
    keys = new  Vector<String>();
    }
    public int GetSize() { return keys.size(); }
    
    
   public int GetId(String key)            // insert if not exist
   {
	 //  System.out.println(key);
	   if(dict.containsKey(key))
		   return dict.get(key);
	   else
	   {
		   int cid = keys.size();
		   keys.add(key);
		   dict.put(key,cid);
		   return cid;
	   }
	}
   
    public int GetIdConst(String key)// return -1 (if not exist)
    {
 	   if(dict.containsKey(key))
		   return dict.get(key);
 	   return -1;
    }


    public String GetKeyWithId( int id) // return "" (if not exist)
    {
        if (id < 0 || id >= keys.size())
            return "";
        return keys.get(id);
    }

    void SaveMappingDict(String file) throws IOException
    {
		BufferedWriter out = new BufferedWriter(new FileWriter( new File(file)));
        for (int i = 0; i < keys.size(); i ++)
            out.write(keys.get(i)+Integer.toString(i)+"\n");			
		out.flush();
		out.close();
    }
    
    void LoadMappingDict(String file)
    {
    	//
    }
    
};