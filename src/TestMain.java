

public class TestMain {

	public static void main(String[] args) {
		System.out.println("Hello, world!\n");

		Planet kerbin = new Planet("KERBIN");
		System.out.println(kerbin.sma);
		kerbin.printInfo(0);
	}

}