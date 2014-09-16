import java.awt.*;
import java.awt.color.ColorSpace;
import java.awt.event.*;
import java.awt.geom.AffineTransform;
import java.awt.image.AffineTransformOp;
import java.awt.image.BufferedImage;
import java.awt.image.BufferedImageOp;
import java.awt.image.ColorConvertOp;
import java.io.*;
import java.util.Vector;
import javax.imageio.ImageIO;
import javax.swing.*;

public class GUI extends JPanel implements MouseListener, MouseWheelListener, MouseMotionListener, KeyListener, ActionListener{

	private static String INPUT_DIR, OUTPUT_DIR;
	private static int FRAME_SIZE = 700;
	private static double THETA_DELTA = Math.toRadians(2);
	private static int RECT_SIZE_DELTA = 2, RECT_TRANSL_DELTA = 2;
	private static final long serialVersionUID = 1012189962657839097L;
	
	private JFrame frame;
	private JLabel label1, label2, label3, label4, label5;
	private JComboBox combo;
	private JTextField textField;
	private JList boxList;
	private JCheckBox autoOCRcheckbox;
	private JButton button1, button2;

	private String[] textlabel = {
			"title", "x-axis-label", "y-axis-label", "x-axis-value", "y-axis-value",
			"z-axis-label", "z-axis-value", "data-value", "legend-title", "legend-category", "misc" };
	private Vector<String> fileNames = new Vector<String>();
	private int fileCounter = 0;
	private BufferedImage currImage = null;
	private Rectangle currImageBounds = null;
	private Dimension currImageOrigDim = null;
	private double aspectR;
	private RotatedTextBox currRect = null;
	private Vector<RotatedTextBox> drawnRect = new Vector<RotatedTextBox>();

	public static void main(String[] args) throws IOException {
		if (args.length < 1) {
			System.err.println("Provide image file directory as commandline argument.");
			System.exit(-1);
		}
		INPUT_DIR = args[0];
		if (args.length < 2) {
			OUTPUT_DIR = INPUT_DIR;
		}
		else {
			OUTPUT_DIR = args[1];
		}
		createDirectoryIfNotExists(INPUT_DIR);
		createDirectoryIfNotExists(OUTPUT_DIR);
		
		GUI gui = new GUI();
		gui.init();
	}
	
	private static void createDirectoryIfNotExists(String dir) {
		File f = new File(dir);
		if (!(f.exists() && f.isDirectory())) {
			f.mkdirs();
		}
	}

	public void init() {
		File folder = new File(INPUT_DIR);
		for (File f : folder.listFiles()) {
			if (f.isFile()) {
				fileNames.add(f.getName());
			}
		}
		this.addMouseListener(this);
		this.addMouseMotionListener(this);
		this.addMouseWheelListener(this);

		frame = new JFrame("VisTextTagger");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		createButtons();
		createCategoryComboBox();
		createTextFieldAndBoxList();
		createLabels();
		frame.getContentPane().add(this);

		setLayout(null); 
		frame.pack();
		frame.setSize(1000, 800);
		frame.requestFocusInWindow();
		frame.addKeyListener(this);
		
		loadImage();
		frame.setVisible(true);
	}

	public void mousePressed(MouseEvent e) {
		// Check if outside image region and reject
		if (e.getX() > currImage.getWidth() || e.getY() > currImage.getHeight()) {
			return;
		}

		if (e.getButton() == MouseEvent.BUTTON1) {
			int x = e.getX();
			int y = e.getY();
			currRect = new RotatedTextBox(x, y, 0, 0, 0, "", "");
			frame.requestFocusInWindow();
			repaint();
		}
	}
	public void mouseDragged(MouseEvent e) {
		checkBoundsAndUpdateRect(e);
	}
	public void mouseReleased(MouseEvent e) {
		checkBoundsAndUpdateRect(e);
	}
	public void mouseWheelMoved(MouseWheelEvent e) {
		checkBoundsAndUpdateRect(e);
	}
	public void mouseEntered(MouseEvent e) {
	}
	public void mouseExited(MouseEvent e) {
	}
	public void mouseClicked(MouseEvent e) {
	}
	public void mouseMoved(MouseEvent e) {
	}

	public void keyReleased(KeyEvent e) {
	}
	public void keyPressed(KeyEvent e) {
		checkBoundsAndUpdateRect(e);
	}
	public void keyTyped(KeyEvent e) {
		if (e.getKeyChar() == 'q') cropCurrRect();
	}
	@Override
	public void paintComponent(Graphics g) { 

		super.paintComponent(g);
		
		Graphics2D g2d = (Graphics2D)g;
		g2d.drawImage(currImage, null, null);

		// Draw the rectangles
		for(RotatedTextBox r : drawnRect) {
			g2d.draw(r.getRotatedRectangle());
		}
		if (currRect != null) {
			g2d.setColor(Color.RED);
			g2d.draw(currRect.getRotatedRectangle());
		}

	}

	/* ----------------------------------- Event response and Saving to file -------------------------------- */
	public void actionPerformed(ActionEvent evt) {

		// Enter pressed in text field
		if (evt.getSource() == textField) {
			saveCurrRectToStack(combo.getSelectedItem());
		}
		// Category chosen from combo box
		else if (evt.getSource() == combo) {
			frame.requestFocusInWindow();
		}
		// Next button pressed
		else if (evt.getSource() == button1){
			saveHandler();
		}
		// Prev button pressed
		else if (evt.getSource() == button2){
			rewindToPreviousImage();
		}
	}
	
	void saveHandler() {
		String fileName[] = fileNames.get(fileCounter).split("\\.");
		File file = new File(OUTPUT_DIR, fileName[0]+".txt");
		
		if (drawnRect.size() > 0) { // Only attempt save when we have data
			if (file.exists()) {
				int b = JOptionPane.showConfirmDialog(null, "Output File exists. Overwrite?");
				if (b == JOptionPane.CANCEL_OPTION || b == JOptionPane.CLOSED_OPTION) {
					return;
				}
				else if (b == JOptionPane.NO_OPTION) {
					clearDataAndProceed();
					return;
				}
			}
			// File does not exist or user has no objections to overwriting
			saveAllRectsToOutputFile(file);
		}
		clearDataAndProceed();
	}
	void clearDataAndProceed() {
		// Clear and proceed to next Image
		clearData();
		fileCounter++;
		if (fileCounter >= fileNames.size()){
			System.exit(0);
			System.out.println("Images over");
		}
		loadImage();
	}
	void rewindToPreviousImage() {
		if (fileCounter > 0) {
			// Clear and go to previous Image
			fileCounter--;
			clearData();
			loadImage();
		}
	}
	
	/** Check bounds after InputEvent e and update current rectangle.
	 * Also perform OCR if Auto-OCR enabled.
	 * @param e
	 */
	void checkBoundsAndUpdateRect(InputEvent e) {
		RotatedTextBox checkRect = new RotatedTextBox(currRect);
		
		if (e instanceof MouseWheelEvent) {
			int notches = ((MouseWheelEvent)e).getWheelRotation();
			double theta = currRect.getTheta() + notches*THETA_DELTA;
			theta %= 2*Math.PI;
			checkRect.setRotation(theta);
		}
		else if (e instanceof MouseEvent) {
			MouseEvent m = (MouseEvent)e;
			int x = m.getX();
			int y = m.getY();
			checkRect.setSize(x - currRect.x,
					y - currRect.y);
		}
		else if (e instanceof KeyEvent) {
			int n = ((KeyEvent)e).getKeyCode();
			
			int mult = 1;
			if (e.isShiftDown()) mult = 5;
			if (e.isControlDown()) mult = 20;
			
			if (n == KeyEvent.VK_W)	 			checkRect.grow(0,mult*RECT_SIZE_DELTA);
			else if (n == KeyEvent.VK_S) 		checkRect.grow(0,-mult*RECT_SIZE_DELTA);
			else if (n == KeyEvent.VK_D)		checkRect.grow(mult*RECT_SIZE_DELTA,0);
			else if (n == KeyEvent.VK_A)		checkRect.grow(-mult*RECT_SIZE_DELTA,0);

			if (n == KeyEvent.VK_LEFT) 			checkRect.translate(-mult*RECT_TRANSL_DELTA, 0);
			else if (n == KeyEvent.VK_RIGHT) 	checkRect.translate(mult*RECT_TRANSL_DELTA, 0);
			else if (n == KeyEvent.VK_UP) 		checkRect.translate(0,-mult*RECT_TRANSL_DELTA);
			else if (n == KeyEvent.VK_DOWN) 	checkRect.translate(0,mult*RECT_TRANSL_DELTA);
			
			else if (n == KeyEvent.VK_ENTER || n == KeyEvent.VK_T) {
				saveCurrRectToStack(combo.getSelectedItem());
			}
			else if (n == KeyEvent.VK_G) {
				combo.setSelectedIndex(Math.max(combo.getSelectedIndex()-1, 0));
			}
			else if (n == KeyEvent.VK_B) {
				combo.setSelectedIndex(Math.min(combo.getSelectedIndex()+1,combo.getItemCount()-1));
			}
			else if (n == KeyEvent.VK_Z) {
				saveHandler();
			}
		}

		checkRect.ensurePositiveSizeAndUpdateRotation();

		if (currImageBounds.contains(checkRect.getBounds())) {
			currRect = checkRect;
		}
		
		if (autoOCRcheckbox.isSelected()) {
			cropCurrRect();
		}
		
		repaint();
		frame.requestFocusInWindow();
	}
	/** OCR text within passed image and set to text field.
	 * @param image 
	 */
	public void ocrExtraction(BufferedImage image) {
		File outputfile = new File("temp.png");
		outputfile.deleteOnExit();
		String ocrText = "", currLine = "";
		try {
			ImageIO.write(image, "png", outputfile);

			// System call to Tesseract OCR
			Runtime r = Runtime.getRuntime();
			Process p = r.exec("tesseract temp.png ocrText -psm 6");
//			Process p = r.exec("tesseract temp.png ocrText");
			p.waitFor();

			// Read text file generated by tesseract
			File f = new File("ocrText.txt");
			f.deleteOnExit();
			BufferedReader br = new BufferedReader(new FileReader(f));
			
			while ((currLine = br.readLine()) != null) {
				ocrText += (currLine + " ");
			}
			if(ocrText.trim().isEmpty()) {
				ocrText = "OCR_FAIL";
			}
			br.close();
		}
		catch (Exception e) {
			e.printStackTrace();
		}

		textField.setText(ocrText.trim());
		textField.requestFocus();
		textField.selectAll();
	}
	/**
	 * Crop out the current rectangle and perform OCR on it.
	 */
	private void cropCurrRect() {
		Rectangle cropRec = currRect.getBounds();
		BufferedImage img = currImage.getSubimage(cropRec.x, cropRec.y, cropRec.width, cropRec.height);
		AffineTransform rotAT = AffineTransform.getRotateInstance(-currRect.getTheta(), cropRec.width, cropRec.height);

		BufferedImageOp bio;
		bio = new AffineTransformOp(rotAT, AffineTransformOp.TYPE_BICUBIC);
		BufferedImage rotImg = bio.filter(img, null);
		ColorSpace cs = ColorSpace.getInstance(ColorSpace.CS_GRAY);
		ColorConvertOp op = new ColorConvertOp(cs, null);
		rotImg = op.filter(rotImg, null);
		
		ocrExtraction(rotImg);
	}
	/**
	 * Save current rectangle and attributes to stack for later output.
	 * @param selectedLabel currently selected label category in the dropdown menu
	 */
	private void saveCurrRectToStack(Object selectedLabel) {
		RotatedTextBox b = new RotatedTextBox(currRect, selectedLabel.toString(), textField.getText());
		drawnRect.add(b);
		boxList.setListData(drawnRect);
		frame.requestFocusInWindow();
		repaint();
	}

	/**
	 * Save all rectangle data to current output file. Returns false if file already exists,
	 * @return For Overwrite dialog: -1 Cancel, 0 Yes - Saved, 1 No - Proceeded with no save
	 */
	private void saveAllRectsToOutputFile(File file) {
		try {
			BufferedWriter outBuffer =
				new BufferedWriter(new FileWriter(file));	

			outBuffer.write("xPixels,yPixels");
			outBuffer.newLine();
			outBuffer.write(currImageOrigDim.width+","+currImageOrigDim.height);
			outBuffer.newLine();
			outBuffer.write("x,y,w,h,theta,category,text");
			outBuffer.newLine();
			for(RotatedTextBox box : drawnRect) {
				outBuffer.write(box.getAsOutputString());
				outBuffer.newLine();
			}
			outBuffer.close();
		}
		catch (IOException e) {
			e.printStackTrace();
		}
	}
	/**
	 * Clear all saved state 
	 */
	public void clearData() {
		drawnRect.clear();
		boxList.setListData(drawnRect);
		currRect = null;
		textField.setText("");
	}
	/**
	 * Load current Image file and resize to Frame size
	 */
	public void loadImage() {
		try {
			currImage = ImageIO.read(new File(INPUT_DIR,fileNames.get(fileCounter)));
			currImageOrigDim = new Dimension(currImage.getWidth(), currImage.getHeight());
			label4.setText("Processing:  " + fileNames.get(fileCounter));
		}
		catch (IOException e) {
			e.printStackTrace();
		}
		// Resize to Frame
		double imgHt = currImage.getHeight();
		double imgWt = currImage.getWidth();
		double wRatio = FRAME_SIZE/imgWt;
		double hRatio = FRAME_SIZE/imgHt;
		aspectR = Math.min(wRatio,hRatio);
		currImage = getScaledInstance(currImage, (int)(aspectR*imgWt), (int)(aspectR*imgHt));
		currImageBounds = new Rectangle(0, 0, currImage.getWidth(), currImage.getHeight());
		repaint();
	}
	/* ----------------------------------- Panel components --------------------------------------------- */
    public BufferedImage getScaledInstance(BufferedImage img, int w, int h) {
        final BufferedImage tmp = new BufferedImage(w, h, BufferedImage.TYPE_INT_ARGB);
        final Graphics2D g2 = tmp.createGraphics();
		g2.setRenderingHint(RenderingHints.KEY_INTERPOLATION,RenderingHints.VALUE_INTERPOLATION_BICUBIC);
        g2.drawImage(img, 0, 0, w, h, null);
        g2.dispose();
        return tmp;
    }
	public void createLabels() {
		label1 = new JLabel("Text Labels", JLabel.CENTER);
		this.add(label1);
		label1.setBounds(800, 110, 150, 20);	
		label1.setVisible(true);
		label1.setBackground(Color.white);
		label1.setOpaque(true);

		label2 = new JLabel("Boxes: Bspace deletes", JLabel.CENTER);
		this.add(label2);
		label2.setBounds(800, 200, 150, 20);	
		label2.setVisible(true);
		label2.setBackground(Color.white);
		label2.setOpaque(true);

		label3 = new JLabel("OCR Output: Q Key", JLabel.CENTER);
		this.add(label3);
		label3.setBounds(800, 20, 150, 20);	
		label3.setVisible(true);
		label3.setBackground(Color.white);
		label3.setOpaque(true);

		label4 = new JLabel("File to be clicked, JLabel.CENTER");
		this.add(label4);
		label4.setBounds(5, FRAME_SIZE + 40, 250, 20);	
		label4.setVisible(true);
		label4.setOpaque(true);
		
		label5 = new JLabel("<html>Keyboard Shortcuts:<br>W,A,S,D = Box Sizing<br>Arrow Keys = Box Translation<br>Shift,Ctrl = Speed Modifiers for above<br>G,B = Change Category Label<br>T or Enter = Save box<br>Q = OCR and select output<br>Z = Save and move to next image</html>");
		this.add(label5);
		label5.setBounds(FRAME_SIZE + 40, 500, 250, 200);	
		label5.setVisible(true);
		label5.setOpaque(true);
	}
	public void createCategoryComboBox() {
		combo = new JComboBox(textlabel);	
		combo.setSelectedIndex(1);
		//menu.setBounds(300, 300, 100, 20);
		this.add(combo);
		combo.setBounds(800, 140, 150, 20);	
		combo.setVisible(true);
		combo.setBackground(Color.white);
		combo.setOpaque(true);
		combo.addActionListener(this);
	}
	public void createTextFieldAndBoxList() {
		textField = new JTextField();
		this.add(textField);
		textField.setBounds(800, 50, 150, 20);
		textField.addActionListener(this);

		boxList = new JList(new DefaultListModel());
		boxList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		boxList.setVisibleRowCount(5);
		boxList.setAutoscrolls(true);
		
		KeyListener boxListKeyListener = new KeyAdapter() {
		    public void keyTyped(KeyEvent e) {
			if ((e.getKeyChar() == '\b') && (drawnRect.size() > 0)) {
			    drawnRect.remove(drawnRect.size()-1);
			    boxList.setListData(drawnRect);
			    repaint();
			}
		    }
		};
		boxList.addKeyListener(boxListKeyListener);
		
		JScrollPane listScroller = new JScrollPane(boxList);
		listScroller.setPreferredSize(new Dimension(250, 200));
		listScroller.setBounds(720, 230, 240, 230);

		this.add(listScroller);
		
		autoOCRcheckbox = new JCheckBox("Auto-OCR");
		autoOCRcheckbox.setBounds(710, 30, 80, 30);
		autoOCRcheckbox.setFocusable(false);
		this.add(autoOCRcheckbox);

	}
	public void createButtons() {
		button1 = new JButton("Next");
		button1.setBounds(870, 700, 100, 30);
		this.add(button1);
		button1.setVisible(true);
		button1.setOpaque(true);
		button1.setBackground(Color.gray);
		button1.addActionListener(this);

		button2 = new JButton("Prev");
		button2.setBounds(750, 700, 100, 30);
		this.add(button2);
		button2.setVisible(true);
		button2.setOpaque(true);
		button2.setBackground(Color.gray);
		button2.addActionListener(this);
	}

	/**
	 * Class encapsulating rotated box with OCR'd text and category label
	 */
	protected class RotatedTextBox extends Rectangle {

		private static final long serialVersionUID = -8980100203661122130L;
		private double theta = 0;
		private String category = "", text = "";
		
		private AffineTransform thetaAT = AffineTransform.getRotateInstance(0);
		
		public RotatedTextBox(int myX, int myY, int myW, int myH, double myTheta, String myCategory, String myText) {
			super(myX, myY, myW, myH);
			setRotation(myTheta);
			setCategory(myCategory);
			setText(myText);
		}
		public RotatedTextBox(RotatedTextBox r) {
			this(r.x, r.y, r.width, r.height, r.theta, "", "");
		}
		public RotatedTextBox(RotatedTextBox r, String myCategory, String myText){
			this(r.x, r.y, r.width, r.height, r.theta, myCategory, myText);
		}
		public RotatedTextBox(Rectangle r, double myTheta) {
			this(r.x, r.y, r.width, r.height, myTheta, "", "");
		}
		
		private void setRotation(double t) {
			theta = t;
			thetaAT.setToRotation(t, getCenterX(), getCenterY());
		}
		
		public void ensurePositiveSizeAndUpdateRotation() {
			if (width <= 0) width=1;
			if (height <= 0) height=1;
			setRotation(theta);
		}
		
		public void draw(Graphics2D g2d) {
			AffineTransform tempAT = g2d.getTransform();
			g2d.setTransform(thetaAT);
			g2d.draw(this);
			g2d.setTransform(tempAT);
		}
		
		public double getTheta() {
			return theta;
		}
		
		public AffineTransform getThetaAT() {
			return new AffineTransform(thetaAT);
		}
		
		@Override
		public Rectangle getBounds() {
			return thetaAT.createTransformedShape(this).getBounds();
		}
		public Shape getRotatedRectangle() {
			return thetaAT.createTransformedShape(this);
		}
		
		public void setCategory(String category) {
			this.category = category;
		}
		public String getCategory() {
			return category;
		}
		public void setText(String text) {
			this.text = text;
		}
		public String getText() {
			return text;
		}
		public String toString() {
			return (text + " - " + category);
		}
		public String getAsOutputString() {
			double nx = (double)(x)/currImage.getWidth();
			double ny = (double)(y)/currImage.getHeight();
			double nw = (double)(width)/currImage.getWidth();
			double nh = (double)(height)/currImage.getHeight();
			return String.format("%.5f,%5f,%5f,%5f,%5f,\"%s\",\"%s\"", 
					nx, ny, nw, nh, theta, category, text);
		}
	}
}