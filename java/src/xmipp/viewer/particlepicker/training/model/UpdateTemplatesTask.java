package xmipp.viewer.particlepicker.training.model;

import javax.swing.SwingWorker;

import xmipp.viewer.particlepicker.training.gui.TemplatesJDialog;


public class UpdateTemplatesTask extends SwingWorker<String, Object>
{

	private static TemplatesJDialog dialog;
	private SupervisedParticlePicker picker;
	private int num;

	public UpdateTemplatesTask(SupervisedParticlePicker picker, int num)
	{
		this.picker = picker;
		this.num = num;
	}
	
	public static void setTemplatesDialog(TemplatesJDialog d)
	{
		dialog = d;
	}

	

	@Override
	protected String doInBackground() throws Exception
	{
		try
		{
			picker.updateTemplates(num);
			
		}
		catch (Exception e)
		{
			throw new IllegalArgumentException(e);
		}
		return "";
	}
	
	 @Override
     protected void done() {
		 if(dialog != null && dialog.isVisible())
				dialog.loadTemplates(true);
     }

}
