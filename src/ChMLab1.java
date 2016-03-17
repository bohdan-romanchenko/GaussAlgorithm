import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

import static java.lang.System.out;

/**
 * Created by Romanchenko Bohdan on 10.03.16.
 */

public class ChMLab1 {

    private static int sizeOfInputArray;

    public static void main(String args[]){
        LSys<Float, Eq> inputArrForEq = null;
        try {
            inputArrForEq = readFromFile(args[0]);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        float[][] inputtedMatrix = printSystem(inputArrForEq);
        int i, j;
        LabAlg<Float, Eq> algorithm = new LabAlg<>(inputArrForEq);
        try{
            algorithm.count();
        }catch (NullPointerException | ArithmeticException e){
            e.printStackTrace();
        }
        Float [] outputVectorArray = new Float[sizeOfInputArray];
        assert inputArrForEq != null;
        for(i = inputArrForEq.size() - 1; i >= 0; i--) {
            Float summary = 0.0f;
            for(j = inputArrForEq.size() - 1; j > i; j--)
                summary = summary + inputArrForEq.valueOf(i, j) * outputVectorArray[j];
            outputVectorArray[i] = (inputArrForEq.valueOf(i, inputArrForEq.size()) - summary) / inputArrForEq.valueOf(i, j);
        }
        printSystem(inputArrForEq);
        printVector(outputVectorArray);
        out.printf("Веткор нев’язків : %s%n", Arrays.toString(getResidual(outputVectorArray, inputtedMatrix)));
    }

    private static Float[] getResidual(Float[] x, float[][] inArray){
        Float[] residual = new Float[sizeOfInputArray];
        Float[] AX = new Float[sizeOfInputArray];
        for (int i = 0; i < sizeOfInputArray; i++){
            AX[i] = 0f;
            for (int j = 0; j < sizeOfInputArray; j++)
                AX[i] = AX[i] + inArray[i][j] * x[j];
            residual[i] = inArray[i][sizeOfInputArray] - AX[i];
        }
        out.println();
        return residual;
    }

    private static LSys<Float, Eq> readFromFile(String path) throws FileNotFoundException {
        LSys<Float, Eq> inputtedList = new LSys<>();
        int i;
        try (Scanner read = new Scanner(new File(path))) {
            sizeOfInputArray = read.nextInt();
            for (i = 0; i < sizeOfInputArray; i++) {
                Eq eq = new Eq();
                eq.readFromFileEq(read);
                inputtedList.inPush(eq);
            }
        }
        return inputtedList;
    }

    private static float[][] printSystem(LSys<Float, Eq> eq){
        float [][] returnArray = new float[sizeOfInputArray][sizeOfInputArray + 1];
        for (int i = 0; i < eq.size(); i++){
            Eq temp = eq.take(i);
            String s = "";
            for (int j = 0; j < temp.size(); j++){
                s = s + String.format("%f; %s", eq.valueOf(i, j), "\t");
                returnArray[i][j] = eq.valueOf(i, j);
            }
            out.println(s);
        }
        out.println();
        return returnArray;
    }

    private static void printVector(Float[] x){
        String s = "";
        for (int i = 0; i < x.length; i++)
            s = s + String.format("x%d = %f; ", i + 1, x[i]);
        out.println(s);
    }

    public interface Gauss<N extends Number, T extends Gauss<N, T>> {
        void inEq(T item);
        void mul(N coefficient);
        N countCoef(N a, N b);
        N in(int index);
        int size();
    }


    public static class LSys<N extends Number, T extends Gauss<N, T>> {
        protected final List<T> list = new ArrayList<>();

        protected T take(int index){
            return list.get(index);
        }

        protected void inPush(T elem){
            list.add(elem);
        }

        protected int size(){
            return list.size();
        }

        protected N valueOf(int i, int j){
            return list.get(i).in(j);
        }
    }

    private static class Eq implements Gauss<Float, Eq> {
        protected final List<Float> equ = new ArrayList<>();
        protected List<Float> getEqu(){
            return equ;
        }
        protected void readFromFileEq(Scanner read) {
            this.equ.clear();
            for (int j = 0; j < sizeOfInputArray + 1; j++)
                if (read.hasNext()) {
                    this.equ.add(read.nextFloat());
                }
        }
        @Override
        public int size(){
            return equ.size();
        }
        @Override
        public void inEq(Eq item){
            ListIterator<Float> i = equ.listIterator();
            ListIterator<Float> j = item.getEqu().listIterator();
            for(; i.hasNext() && j.hasNext();){
                Float a = i.next();
                Float b = j.next();
                i.set(a + b);
            }
        }
        @Override
        public void mul(Float coefficient){
            for(ListIterator<Float> i = equ.listIterator(); i.hasNext();){
                Float futureVal = i.next();
                i.set(futureVal * coefficient);
            }
        }
        @Override
        public Float countCoef(Float valA, Float valB){
            if (valA == 0.0f) return 1.0f;
            return -valB/valA;
        }
        @Override
        public Float in(int index){
            return equ.get(index);
        }
    }

    protected static class LabAlg<Number extends java.lang.Number, eq extends Gauss<Number, eq>> {
        protected final ThreadLocal<LSys<Number, eq>> list = new ThreadLocal<LSys<Number, eq>>() {
            @Override
            protected LSys<Number, eq> initialValue() {
                return null;
            }
        };
        protected LabAlg(LSys<Number, eq> lSys){
            list.set(lSys);
        }

        protected void count() throws ArithmeticException {
            if (!checkSystem(list.get())) {
                throw new ArithmeticException("Incorrect system for Gauss Method");
            }
            for(int i = 0; i < list.get().size() - 1; i++)
                for(int j = i + 1; j < list.get().size(); j++){
                    Number k = list.get().take(j).countCoef(list.get().take(j).in(i), list.get().take(i).in(i));
                    list.get().take(j).mul(k);
                    list.get().take(j).inEq(list.get().take(i));
                }
        }

        protected boolean checkSystem(LSys<Number, eq> system){
            if (2 > system.size()) return false;
            for(int i = 0; i < system.size(); i++)
                if (system.take(i).size() != (system.size() + 1))
                    return false;
            return true;
        }
    }
}
