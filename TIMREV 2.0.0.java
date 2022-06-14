/* TIMREV 2.0.0 */
import java.util.*;
import java.lang.Math.*;
class Array_Initial_Tools {
    public double[] Array_Reset(double[] A, int length) { /*Resets given Array*/
        for (int i = 0; i < length; i++) {
            A[i] = 0;
        }
        return A;
    }
    public double [][]Matrix_Reset(double [][]M,int rows,int columns){ /*Resets given Matrix*/
        for(int i=0;i<rows;i++){
            for(int j=0;j<columns;j++){
                M[i][j] = 0;
            }
        }
        return M;
    }
    public double[] Select_Rows_From_Matrix(double[][] M, int columns, int n) { /*Returns an Array from a given Matrix*/
        double[] R = new double[columns];
        Array_Reset(R, columns);
        for (int i = 0; i < columns; i++) {
            R[i] = M[n][i];
        }
        return R;
    }
    public double [][] Set_Matrix_Array_Trace (double [][]M,int rows, int columns){
        Matrix_Reset(M,rows,columns);               
        for(int i=0;i<rows;i++){
            M[i][0]=(i+1);
        }
        return M;
    }
}
class Member{
    public double [] Member_Initial_Data (double t, double d, double L, double K){ /*Gets Initial Member Data*/
        double [] MID = new double[7];
        Array_Initial_Tools AIT = new Array_Initial_Tools();
        AIT.Array_Reset(MID,7);
        MID[0] = t * d; /*Cross section Area */
        MID[1] = (t *(Math.pow(d,3)))/12;  /*Moment of Inertia around x Axis*/
        MID[2] = (d * (Math.pow(t,3)))/12; /*Moment of Inertia around y Axis*/
        MID[3] = (Math.sqrt(MID[1]/MID[0])); /*The gyration radius around x Axis*/
        MID[4] = (Math.sqrt(MID[2]/MID[0])); /*The gyration radius around y Axis*/
        MID[5] = (K * L * 100)/d; /*Slenderness ratio around x Axis*/
        MID[6] = (K * L * 100)/t; /*Slenderness ratio around y Axis*/
        return MID;
    }
    public double Member_Linear_Capacity (double t,double d,double L,double K,double n){/*Returns Various Types of Member Capacity in Linear(Elastic) State*/
        double [][] M = new double[][]
        double [] R = new double[6];
        double [] MID = new double[7];
        Array_Initial_Tools AIT = new Array_Initial_Tools();
        AIT.Array_Reset(MID,7);
        MID = Member_Initial_Data(t,d,L,K);
    }
}
public class TIMREV {
    public static void main (String []args){
        Scanner scan = new Scanner(System.in);
        System.out.println("Hello world");
    }
}
