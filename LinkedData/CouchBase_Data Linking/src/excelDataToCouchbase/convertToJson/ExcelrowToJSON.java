package excelDataToCouchbase.convertToJson;

import java.text.NumberFormat;
import java.text.ParsePosition;
import java.util.ArrayList;
import java.util.List;

import com.couchbase.client.java.document.JsonDocument;
import com.couchbase.client.java.document.json.JsonObject;
import com.incesoft.tools.excel.xlsx.Cell;
import com.incesoft.tools.excel.xlsx.Sheet;
import com.incesoft.tools.excel.xlsx.Sheet.SheetRowReader;

import excelDataToCouchbase.connect.ConnectUpload;


public class ExcelrowToJSON {
	
	ConnectUpload upload = new ConnectUpload();
	static List<Object> spaltenName;
	
	public void convert(Sheet sheet, String filename) throws Exception{
		List<JsonDocument> list = new ArrayList<JsonDocument>();
		SheetRowReader reader = sheet.newReader();
		Cell[] row;
		int rowPos = 0;
		row = reader.readRow();
		
		spaltenName = new ArrayList<Object>();
		
		for (Cell cell:row){
			spaltenName.add(cell.getValue());
		}
		rowPos++;
		JsonObject content = JsonObject.empty();
		while((row = reader.readRow())!=null){
			int i=0;
			for(Cell cell:row){
				if(cell == null){
					content.put((String) spaltenName.get(i), "");
				}else if(isNumeric(cell.getValue())){
					if(!cell.getValue().equals("")){
						if(cell.getValue().contains(".")){
							content.put((String) spaltenName.get(i), Double.parseDouble(cell.getValue().toString()));
						} else{
							content.put((String) spaltenName.get(i), NumberFormat.getInstance().parse(cell.getValue()));
						}
						
					} else{
						content.put((String) spaltenName.get(i), "");
					}
				} else{
					content.put((String) spaltenName.get(i), cell.getValue());
				}
				i++;
			}
			JsonDocument doc = JsonDocument.create(filename + "_" + rowPos, content);
			list.add(doc);
//			upload.connectAndUpload(doc);
			rowPos++;
			if(rowPos==2){
				System.out.print("IN PROGRESS: .");
			}
			if(rowPos%100 == 0)
			System.out.print(".");
		}
		upload.connectAndUpload(list);
		upload.close();

	}
	
	public static boolean isNumeric(String str){
	  NumberFormat formatter = NumberFormat.getInstance();
	  ParsePosition pos = new ParsePosition(0);
	  formatter.parse(str, pos);
	  return str.length() == pos.getIndex();
	}
}
