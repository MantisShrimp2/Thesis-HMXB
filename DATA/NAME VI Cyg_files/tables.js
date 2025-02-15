
$(document).ready(function() {
//	$('table.sortable').each(function() {
//		$(this).alternateRowColors();
//	});

//	$('#tab1').find('tr').eq(2).find('td').addClass('highlight');
//	$('#tab1').find('tr').eq(2).find('td').eq(1).addClass('highlight');
//	$('#tab2').find('tr').find('td:even').addClass('highlight');

	$('table.sortable').each(function() {
		$(this).sortIt();
	});
});
