import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.lang.Math;

public class Planet {

/*
	The Planet class 
*/

	// Class variables
	private static final String FILENAME = "solar_bodies.csv";
	public String name, parent;
	public double sma, inc, ecc, lpe, lan, mna;
	public double period, meanmotion, mu, radius, soi;
	public double orb_energy, semiparam;

/* 	VARIABLE DEFINITIONS
	name 		Name of body, lowercase string
	parent		Parent of body, lowercase string

	sma			Semimajor axis, kilometers
	inc 		Inclination, degrees
	ecc			Eccentricity, unitless
	lpe			Argument of periapsis, degrees
	lan			Longitude of ascending node (actually right ascension), degrees
	mna 		Mean anomaly at 0.0 UT, radians

	period		Sidereal orbital period, seconds
	meanmotion 	Mean orbital motion, radians per second

	mu 			Standard gravitational parameter, m^3/s^2
	radius		Equatorial radius, meters
	soi 		Sphere of influence, meters

	orb_energy 	Orbital energy, m^2/s^2
	semiparam 	Semilatus rectum, m
*/

	// Constructor
	public Planet(String name) {

		// Sanitize input and initialize readers
		name = name.toLowerCase();
		BufferedReader br = null;
		FileReader fr = null;

		// Wrap full file read in a try/catch block
		try {

			fr = new FileReader(FILENAME);
			br = new BufferedReader(fr);

			String current_line;

			// Read the file. Check for 
			while ((current_line = br.readLine()) != null) {
				// Strip whitespace and non-visible characters
				current_line.replaceAll("\\s+","");
				// Split string assuming comma delimited
				String[] parts = current_line.split(",");

				// Check if first element of line matches user-given body
				if (name.equals(parts[0].toLowerCase())) {

					this.name = name;
					this.parent = parts[1].toLowerCase();

					this.sma = Double.parseDouble(parts[2]);
					this.inc = Double.parseDouble(parts[3]);
					this.ecc = Double.parseDouble(parts[4]);
					this.lpe = Double.parseDouble(parts[5]);
					this.lan = Double.parseDouble(parts[6]);
					this.mna = Double.parseDouble(parts[7]);

					this.period = Double.parseDouble(parts[8]);
					this.meanmotion = Double.parseDouble(parts[9]);

					this.mu = Double.parseDouble(parts[10]);
					this.radius = Double.parseDouble(parts[11]);

					this.orb_energy = -1*this.mu / (2*this.sma);
					this.semiparam = this.sma * (1 - pow(this.ecc,2));

					// Planet found, stop reading
					break;
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			// Close file and buffer readers
			try {
				if (br != null)
					br.close();
				if (fr != null)
					fr.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}

		}

	}


	// Return Keplerian elements at time t
	public double[] getKepler(double t) {
		double kepler_vec[] = new double[6];

		kepler_vec[0] = this.sma;
		kepler_vec[1] = this.inc;
		kepler_vec[2] = this.ecc;
		kepler_vec[3] = this.lpe;
		kepler_vec[4] = this.lan;
		kepler_vec[5] = getTrueAnomaly(this.mna + t*this.meanmotion, this.sma, this.ecc);

		return kepler_vec;
	}


	// Compute true anomaly from mean anomaly
	public static double getTrueAnomaly(double mean_anom, double a, double e) {
		// Find coterminal mean anomaly between 0 and 2*pi
		if (mean_anom > 2*Math.PI) {
			mean_anom = getCoterm(mean_anom);
		}

		// Use Newton-Raphson iteration to compute eccentric anomaly
		// Use mean anomaly as initial guess (true for circular)
		double ecc_anom = mean_anom;
		double ecc_prev;	// Previous eccentric anomaly for iterations
		double r;	// Magnitude of position vector
		double sin_true;	// Y component of position vector
		double cos_true;	// X componenet of position vector
		double result;

		// Set error tolerance
		double tol = 0.0000001;

		// Set error to large number
		double err = tol + 1;

		// If circular, true anomaly is equal to mean anomaly
		if (e == 0) return mean_anom;

		// Iterate eccentric anomaly
		while (err >= tol) {
			ecc_prev = ecc_anom;
			ecc_anom = ecc_anom + (mean_anom - ecc_anom + e*Math.sin(ecc_anom)) / (1 - e*Math.cos(ecc_anom));
			err = Math.abs(ecc_anom-ecc_prev);
		}

		// Magnitude of position vector
		r = a*(1 - e*Math.cos(ecc_anom));
		// Components of position vector
		sin_true = (Math.sin(ecc_anom)*Math.sqrt(1-Math.pow(e,2))) / (1-e*Math.cos(ecc_anom));
		cos_true = (Math.cos(ecc_anom)-e) / (1-e*Math.cos(ecc_anom));

		result = Math.atan2(sin_true,cos_true);
		if (result < 0) result += 2*Math.PI;

		return result;
	}


	// Get inertial position and velocity in Cartesian coordinates at time t
	public double[] getState(double t) {
		// Local variables
		double state[6];
		double mean_anom, true_anom;
		double energy, semiparam;

		// Compute true anomaly
		mean_anom = this.mna + t*this.meanmotion;
		true_anom = getTrueAnomaly(mean_anom, this.sma, this.ecc);

		// Check for edge cases
		if (this.ecc == 0) this.lpe = 0;	// Circular
		if (this.inc == 0) this.lan = 0; 	// Equatorial

		// Construct coordinates in PQW frame
		double[] Rpqw = {	this.semiparam*Math.cos(true_anom)/(1+this.ecc*Math.cos(true_anom)),
							this.}
	
		return state;
	}


	// Print variables at time t to the console
	public void printInfo(double t) {
		// Local variables
		double mean_anom, true_anom;

		// Compute mean and true anomaly at time t
		mean_anom = this.mna + t*this.meanmotion;
		true_anom = getTrueAnomaly(mean_anom, this.sma, this.ecc);

		// Print constant orbital information
		System.out.println("SMA (m), ECC, INC (deg), LPE (deg), LAN (deg)");
	}


	// Find coterminal angle between 0 and 2*pi
	public static double getCoterm(double theta) {
		// Local variables
		double result = theta;

		while (result > 1000*2*Math.PI) {
			result -= 1000*2*Math.PI;
		}

		while (result > 2*Math.PI) {
			result -= 2*Math.PI;
		}

		return result;
	}


}