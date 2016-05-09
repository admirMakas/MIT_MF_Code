import com.phoenix_int.ModelCenter.*;
public class Main {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		try
		{
		ModelCenter model = new ModelCenter();
		model.loadModel("C:/Users/Cory/Desktop/NASA Optimizer/Rosenbrock.pxc");
		
		for (int i=0; i<(args.length-1); i++){
			String numString = Integer.toString(i);
			//ModelCenter Model must have variables named x0 through x(n-1)
			model.setValue("Model.Matlab.x" + numString, args[i]);
		}
		
		//ModelCenter Model must have fidelity variable named "fidelity"
		model.setValue("Model.Matlab.fidelity", args[args.length-1]);
		model.run(null);
		System.out.println(model.getValue("Model.Matlab.fnEval"));
		}
		catch ( Exception e )
		{
		System.out.println( "Error: " + e.getMessage() );
		}		
	}
}
