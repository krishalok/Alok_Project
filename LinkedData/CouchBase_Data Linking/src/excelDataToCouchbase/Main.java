package excelDataToCouchbase;

import java.awt.Container;
import java.io.File;
import java.util.Scanner;

import com.incesoft.tools.excel.ReaderSupport;
import com.incesoft.tools.excel.support.XLSXReaderSupport;
import com.incesoft.tools.excel.xlsx.Sheet;
import com.incesoft.tools.excel.xlsx.SimpleXLSXWorkbook;
import com.sun.glass.ui.Cursor;

import excelDataToCouchbase.convertToJson.ExcelrowToJSON;


public class Main {
	public static void main(String[] args){
		
		
		System.out.println("Welcome! You´re using a Excel (xlsx) importer for Couchbase!\nPlease give the fullpath of your Excel (.xlsx) File: ");
		
		Scanner scanner = new Scanner(System.in);
		String filepath = scanner.nextLine();
		File myFile = new File(filepath);
		SimpleXLSXWorkbook workbook = new SimpleXLSXWorkbook(myFile);
		Sheet mySheet = workbook.getSheet(0, false);
		
		String[] split = filepath.split("\\\\");
		String filename = split[split.length - 1];
		String foldername = filename.substring(0, filename.length()-5);
		ExcelrowToJSON excelToJson = new ExcelrowToJSON();
		try{
			long startTime = System.nanoTime();
			excelToJson.convert(mySheet, foldername);
			long endTime = System.nanoTime();
			System.out.println();
			System.out.println(filename  + ": upload successful finished!");
			System.out.println("Duration: " + (endTime - startTime));
		}catch (Exception e) {
			e.printStackTrace();
		}
		scanner.close();
	}
}
