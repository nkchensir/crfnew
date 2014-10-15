package crf;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class test {
     private static int a;
	public void tt(int[] a)
	{
	 ta(a);
	}
	public void ta(int[] a)
	{
	 a[0]=10000;
	}
		
	public void ttt(int b)
	{
	 b=10000;
	 a=100;
	}
	public void td(testdf d)
	{
	 d.a=100000000;
	 a=100;
	}
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		String[] a=null;
		a="1 2 3 4 ".split(" ");
		System.out.println(Math.abs(-2.73));		
		test C= new test();
		C.a  = 10;
		boolean b;
		b = true;
		C.ttt(C.a);
	     System.out.println(C.a);		
        ArrayList L = new ArrayList();
        int[] A;
        int [] B = new int[4];
        A = new int[4];
        A[0] = 1;
        A[1] = 2;
        A[2] = 3;
        A[3] = 4;

        test T = new test();
        T.tt(A);
        System.out.println("A0:"+A[0]);    
        testdf dd= new testdf();
        testdfs dfs= new testdfs();
        dfs.b=1234;
        System.out.println(((testdf) dfs).a);
        System.out.println(dfs.b);
        testdf[] ds = new testdf[3];
        ds[0] = new testdf();
        ArrayList L1 = new ArrayList();
        ds[0].a = 31;
        L1.add(ds[0]);
        ((testdf) L1.get(0)).a = 24;
        
        System.out.println(ds[0].a);
        /*
        C.td(dd);
	     System.out.println(dd.a);	      
        L1.add(A);
        System.out.println(L1.size());
        System.out.println(A.length);
        int w = ((int[]) L1.get(0))[0];
        w =12;
        System.out.println(((int[]) L1.get(0))[0]);
        B[2] =1000;
        System.out.println(A[2]);
        L1.indexOf(A);*/
	}

}
