import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * Created by Romanchenko Bohdan on 10.03.16.
 */
public class ChMLab1 {

    private static int DEFAULT_EQUATIONS_NUMBER;

    public static void main(String args[]){
        LinearSystem<Float, MyEquation> list = null;
        try {
            list = readFromFile(args[0]);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        float[][] array = printSystem(list);
        int i, j;
        Algorithm<Float, MyEquation> alg = new Algorithm<>(list);
        try{
            alg.calculate();
        }catch (NullPointerException | ArithmeticException e){
            System.out.println(e.getMessage());
            System.exit(0);
        }
        Float [] x = new Float[DEFAULT_EQUATIONS_NUMBER];
        assert list != null;
        for(i = list.size() - 1; i >= 0; i--) {
            Float sum = 0.0f;
            for(j = list.size() - 1; j > i; j--) {
                sum += list.itemAt(i, j) * x[j];
            }
            x[i] = (list.itemAt(i, list.size()) - sum) / list.itemAt(i, j);
        }
        printSystem(list);
        printVector(x);
        System.out.println("Веткор нев’язків : " + Arrays.toString(getResidual(x, array)));
    }

    public static Float[] getResidual(Float[] x, float[][] inArray){
        Float[] residual = new Float[DEFAULT_EQUATIONS_NUMBER];
        Float[] AX = new Float[DEFAULT_EQUATIONS_NUMBER];
        for (int i = 0; i < DEFAULT_EQUATIONS_NUMBER; i++){
            AX[i] = 0f;
            for (int j = 0; j < DEFAULT_EQUATIONS_NUMBER; j++)
                AX[i] += inArray[i][j] * x[j];
            residual[i] = inArray[i][DEFAULT_EQUATIONS_NUMBER] - AX[i];
        }
        System.out.println();
        return residual;
    }

    public static LinearSystem<Float, MyEquation> readFromFile(String path) throws FileNotFoundException {
        LinearSystem<Float, MyEquation> list = new LinearSystem<>();
        int i;
        Scanner read = new Scanner(new File(path));
        DEFAULT_EQUATIONS_NUMBER = read.nextInt();
        for (i = 0; i < DEFAULT_EQUATIONS_NUMBER; i++){
            MyEquation eq = new MyEquation();
            eq.readFromFileEq(read);
            list.push(eq);
        }
        return list;
    }

    public static float[][] printSystem(LinearSystem<Float, MyEquation> system){
        float [][] returnArray = new float[DEFAULT_EQUATIONS_NUMBER][DEFAULT_EQUATIONS_NUMBER + 1];
        for (int i = 0; i < system.size(); i++){
            MyEquation temp = system.get(i);
            String s = "";
            for (int j = 0; j < temp.size(); j++){
                s += String.format("%f; %s", system.itemAt(i, j), "\t");
                returnArray[i][j] = system.itemAt(i, j);
            }
            System.out.println(s);
        }System.out.println();
        return returnArray;
    }

    public static void printVector(Float [] x){
        String s = "";
        for (int i = 0; i < x.length; i++){
            s += String.format("x%d = %f; ", i + 1, x[i]);
        }System.out.println(s);
    }

    public interface Gauss<N extends Number, T extends Gauss<N, T>> {
        void addEquation(T item);
        void mul(N coefficient);
        N findCoefficient(N a, N b);
        N at(int index);
        int size();
    }


    public static class LinearSystem<N extends Number, T extends Gauss<N, T>> {
        private List<T> list = new ArrayList<>();

        public T get(int index){
            return list.get(index);
        }

        public void push(T elem){
            list.add(elem);
        }

        public int size(){
            return list.size();
        }

        public N itemAt(int i, int j){
            return list.get(i).at(j);
        }
    }

    public static class MyEquation implements Gauss<Float, MyEquation> {
        private List<Float> equation = new ArrayList<>();
        public List<Float> getEquation(){
            return equation;
        }
        public void readFromFileEq(Scanner read) throws FileNotFoundException {
            this.equation.clear();
            for (int j = 0; j < DEFAULT_EQUATIONS_NUMBER + 1; j++)
                if (read.hasNext())
                    this.equation.add(read.nextFloat());
        }
        @Override
        public int size(){
            return equation.size();
        }
        @Override
        public void addEquation(MyEquation item){
            ListIterator<Float> i = equation.listIterator();
            ListIterator<Float> j = item.getEquation().listIterator();
            for(; i.hasNext() && j.hasNext();){
                Float a = i.next();
                Float b = j.next();
                i.set(a + b);
            }
        }
        @Override
        public void mul(Float coefficient){
            for(ListIterator<Float> i = equation.listIterator(); i.hasNext();){
                Float next = i.next();
                i.set(next * coefficient);
            }
        }
        @Override
        public Float findCoefficient(Float a, Float b){
            if (a == 0.0f) return 1.0f;
            return -b/a;
        }
        @Override
        public Float at(int index){
            return equation.get(index);
        }
    }

    public static class Algorithm<N extends Number, T extends Gauss<N, T>> {
        LinearSystem<N, T> list = null;
        public Algorithm(LinearSystem<N, T> system){
            list = system;
        }

        public void calculate() throws NullPointerException, ArithmeticException{
            if (!checkSystem(list)){
                throw new ArithmeticException("Incorrect system for Gauss Method");
            }
            for(int i = 0; i < list.size() - 1; i++){
                for(int j = i + 1; j < list.size(); j++){
                    N k = list.get(j).findCoefficient(list.get(j).at(i), list.get(i).at(i));
                    list.get(j).mul(k);
                    list.get(j).addEquation(list.get(i));
                }
            }
        }

        private boolean checkSystem(LinearSystem<N, T> system){
            if (system.size() < 2) return false;
            for(int i = 0; i < system.size(); i++){
                if (system.get(i).size() != (system.size() + 1)){
                    return false;
                }
            }
            return true;
        }
    }
}
