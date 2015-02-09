/**
 * @author George Dietz
 * 
 * That which makes it go
 */

$(document).ready(function() {
	$( "#from" ).datepicker({
		defaultDate: "-1m",
		changeMonth: true,
		changeYear: true,
		numberOfMonths: 1,
		onClose: function( selectedDate ) {
			$( "#to" ).datepicker( "option", "minDate", selectedDate );
		}
	}).datepicker("setDate", "-1m");

	$( "#to" ).datepicker({
		defaultDate: "+0",
		changeMonth: true,
		changeYear: true,
		numberOfMonths: 1,
		onClose: function( selectedDate ) {
			$( "#from" ).datepicker( "option", "maxDate", selectedDate );
		}
	}).datepicker("setDate", "+0");
	

	
	var main_table = $('#main_table').dataTable( {
		"aaSorting": [],
		"aoColumns": [
		{ "sWidth": "11em"},
		null,
		null,
		null,
		],
		"iDisplayLength": 10,
	} );
	go_handler();

	/* start of btn handling function */
	$( "#go_btn" ).click(go_handler); /* end of btn handling function */
	
	
	function go_handler() {
		var from_date = $('#from').datepicker("getDate");
		var to_date = $('#to').datepicker("getDate");
		
		from_date = $.datepicker.formatDate("yymmdd", from_date);
		to_date = $.datepicker.formatDate("yymmdd", to_date);
		
		var urltoload = '/'+from_date+'/'+to_date+'/';
		
		var json_data;
		$('body').addClass("loading");
		$.getJSON(urltoload,function(json){
			json_data = json;
			main_table.fnClearTable();
			main_table.fnAddData(json_data.aaData);

			$('body').removeClass("loading");
			main_table.fnAdjustColumnSizing();              
		});
	}
	
	
	$('#export_csv').click(function (event) {
		var csv = get_csv(main_table);
		var csvData = 'data:application/csv;charset=utf-8,' + encodeURIComponent(csv);

		$(this)
			.attr({
			'download': 'pmids.csv',
				'href': csvData,
				'target': '_blank'
		});

	});
	$('#export_txt').click(function (event) {
		var csv = get_csv(main_table);
		var csvData = 'data:application/plain;charset=utf-8,' + encodeURIComponent(csv);

		$(this)
			.attr({
			'download': 'pmids.txt',
				'href': csvData,
				'target': '_blank'
		});
	});
	
	function get_csv(table) {
		var data = table.fnGetData();

		var pmids = data.map(function(row) {
			return [$.parseHTML(row[1])[0].innerHTML];
		});
		var csv = Papa.unparse(pmids);

		
		return csv;
	}
	
	
	
	
	
	
}); // end of $(document).ready()