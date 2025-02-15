var mytable=null;

jQuery.fn.print = function(message) {
	return this.each(function() {
		$('<div class="result" />')
			.text(String(message))
			.appendTo($(this).find('.result'));
	});
};

jQuery.fn.alternateRowColors = function() {
	$('tbody tr:odd', this).removeClass('even').addClass('odd');
	$('tbody tr:even', this).removeClass('odd').addClass('even');
	return this;
};
jQuery.fn.sortIt = function() {
	var table = $(this);
	$('th', table).each(function(column) {
		var header = $(this);
		var findSortKey;
		if (header.is('.sort-alpha')) {
			findSortKey = function(cell, order) {
				var str = jQuery.trim(cell.text().toUpperCase());
				if (str.indexOf('~') >= 0 || str == '') str = order == 1 ? '~':' ';
				return str;
			};
		} else if (header.is('.sort-int')) {
			findSortKey = function(cell, order) {
				var num;
				var str = jQuery.trim(cell.text());
				if (str.indexOf('~') >= 0 || str == '' || str == ' ' ) num = order*2147483647;
				else num = parseInt(cell.text());
				return num;
			};
		} else if (header.is('.sort-float')) {
			findSortKey = function(cell, order) {
				var str = jQuery.trim(cell.text());
				if (str.indexOf('~') >= 0 || str == '' || str == ' ' ) return order*1.0e+30;
				else return parseFloat(cell.text());
			};
		} else if (header.is('.sort-sexa')) {
			findSortKey = function(cell, order) {
				var str = jQuery.trim(cell.text());
				if (str.indexOf('~') >= 0 || str == '') return order*1.0e+30;
				var digits = str.match(/([0-9.+-]+)/g);
				var num = 0.0;
				var maxidx;
				var sign = 1;
				$.each(digits, function(idx, val) {
					if (num == 0.0 && val.charAt(0) == '-') {
						sign = -1;
						val = val.substring(1,val.length);
					     }
					var f = parseFloat(val);
					num = num*60+f;
					maxidx = idx;
				});
				if (maxidx == 1) { num *= 60.0; }
				if (maxidx == 0) { num *= 3600.0; }
				return sign == 1 ? num : -num;
			};
		} else if (header.is('.sort-sp')) {
			findSortKey = function(cell, order) {
				var str = jQuery.trim(cell.text());
				if (str.indexOf('~') >= 0 || str == '') return order == 1 ? '99~' : '   ';
				var code = 10+'OBAFGKMSCDR'.indexOf(str.substring(0,1));
				if (code < 10) code = 98;
				var key = code+str;
				return key;
			};
                  } else if (header.is('.sort-img')) {
                           findSortKey = function(cell, order) {
                                    var strImg = jQuery.trim(cell.html());
				    var strDjin = jQuery.trim($(cell.context).find("td").eq(3).text());

                                    //if (strImg.indexOf('img')>=0 && strDjin.) return 1;
                                    if (strDjin.indexOf('D')>=0) return 2;
				    else if (strImg.indexOf('img')>=0) return 1;
                                else return 0;
                           };

		  } else if (header.is('.sort-tak')) {
			   findSortKey = function(cell, order) {
				    var str = jQuery.trim(cell.text());
				    if (str.indexOf('T')>=0||str.indexOf('A')>=0||str.indexOf('K') >= 0) 
					return 1;
					else return 0;
			   };
		  } else if (header.is('.sort-djin')) {
		// ((100*$title+30*$abstract+50*$keyword)*$pages^0.25+100*$subtitle+100*$caption+20*$data+100*$figure+50*($nbocc-$title-$abstract-$keyword-$subtitle-$caption-$data-$figure)*($text-$calibrator)/$pages)/$age^0.2
			   findSortKey = function(cell, order) {
				    var str = jQuery.trim($(cell.context).find("td").eq(3).text()+$(cell.context).find("td").eq(2).text());
				var title = str.indexOf('T')>=0?1:0;
				var abst = str.indexOf('A')>=0?1:0;
				var keyword = str.indexOf('K')>=0?1:0;
				var subtitle = str.indexOf('S')>=0?1:0;
				var caption = str.indexOf('C')>=0?1:0;
				var data = str.indexOf('D')>=0?1:0;
				var figure = str.indexOf('F')>=0?1:0;
				var nbocc = parseInt($(cell.context).find("td").eq(1).text());
				if (nbocc>0) nbocc = parseInt(str.substring(str.lastIndexOf(' ')+1,str.length));
				else nbocc = 1; // definit à 1 le nb min d'occurence
				var text = str.indexOf('X')>=0?1:0;
				var calibrator = str.indexOf('*')>=0?1:0;
				var age = 10;
			//	var age = 2014-parseInt(($(cell.context).find("td").eq(1).text()).subtring(0,4));
				var sum = ((100*title+30*abst+50*keyword)+100*subtitle+100*caption+20*data+100*figure+50*(nbocc-title-abst-keyword-subtitle-caption-data-figure)*(text-calibrator))/age^0.2;
				return sum;
			   };

		} else if (header.is('.sort-pm')) {
			findSortKey = function(cell, order) {
				var pm = cell.text().match(/([0-9.~+-]+)/g);
				if (pm[0] == '~' || pm[1] == '~') return order*9.0E+20;
				num = pm[0]*pm[0]+pm[1]*pm[1];
				return num;
			};
		} else if (header.is('.sort-dim')) {
			findSortKey = function(cell, order) {
				var pos = cell.text().indexOf(' ');
				if (pos < 0) return order*9.0E+20;
				var majaxis = jQuery.trim(cell.text().substring(0, pos));
				if (majaxis == '~') return order*9.0E+20;
				return majaxis;
			};
		};
				

		if (findSortKey) {
			header.addClass('clickable').hover(
				function() {header.addClass('hover');},
				function() {header.removeClass('hover');});
			var order;
			header.click(function() {
				var tbody = table.find('tbody');
				var rows = tbody.find('tr').get();	// rows = array
				if (header.is('.order-asc')) {
					order = -1;
					header.removeClass('order-asc').addClass('order-desc');
				} else if (header.is('.order-desc')) {
					order = 1;
					header.removeClass('order-desc').addClass('order-asc');
				} else {
					if (header.is('.order-2-desc')) {
						order = -1;
					}
					else {
						order = 1;
						header.addClass('order-asc');
					}
				}
				$.each(rows, function(idx, row) {
					var cell = $(row).children('td').eq(column);
					row.sortKey = findSortKey(cell, order);	// ???? row.sortKey ?
					row.position = idx;
				});
				rows.sort(function(a,b) {
					if (a.sortKey < b.sortKey) return -order;
					if (a.sortKey > b.sortKey) return order;
					return a.position < b.position ? -1 : 1;
				});
				$.each(rows, function(idx, row) {
					table.children('tbody').append(row);
					row.sortKey = null;
				});
				var thead = table.find('thead');
				var cols = thead.find('th').get();

				$.each(cols, function(idx, col) {
					if (idx != column) {
						$(col).removeClass('order-asc').removeClass('order-desc').addClass('order-gen');
					}
				});
				/* table.alternateRowColors(); */
			});
		}
	});
	return table;
}


var tab2map =function (tab) {
	var res={};
	tab.forEach(function(element) {
		res[element[0]]=element[1];
	});
	return res;
}

// fonction pour eviter de lancer trop souvent une méthode
// utilisée quand on veut faire un appel à chaque touche
var debounce = function(fn, delay) {
	var timer = null;
	return function () {
		var context = this, args = arguments;
		clearTimeout(timer);
		timer = setTimeout(function () {
			fn.apply(context, args);
		}, delay);
	};
};

var getBibcodeTitle = function(elemnt) {
	var cache = sessionStorage.bibcode_cache;
	if (!cache) cache="";
	bib=$(elemnt).text().trim();
	if (! $(elemnt).hasClass('bibcode') || bib.length < 18){
		return;
	}
	if (cache.indexOf(bib)==-1) {
	$.ajax({
		url: 'https://simbad.u-strasbg.fr/simbad/sim-tap/sync',
		type:'GET',
		cache:true,
		data: {query: 'select bibcode,title from ref where bibcode=\''+bib+'\'', format: 'text', lang: 'adql', request :'doQuery'},
		async:false
	})
	// une fois qu'on a le resultat, remplit le titre de l'element
	.done(function(result) {
		var b = result.split("\n")[2].split("|")[0].replace(/"/g,'');
		var n = result.split("\n")[2].split("|")[1];
		if (n.length > 2) {
			$(elemnt).attr('title',n);
			// fabrique un nouveau cache avec le nouveau titre
			var newcache = {};
			if (cache.length>0)
				newcache = JSON.parse(cache);
	
			newcache[b]=n;
			sessionStorage.bibcode_cache=JSON.stringify(newcache);
		}
		else {
			console.log("erreur : reponse tap Null : "+result);
		}

	});
	}
	else {
		var cacheBib = JSON.parse(cache);
		
		$(elemnt).attr('title',cacheBib[bib]);
	}
	

};
var getBibcodeTitles = function(bibs, cachebib) {
	var newcache = {};
	$.ajax({
		url: '/simbad/sim-tap/sync',
		type:'GET',
		cache:true,
		data: {query: 'select bibcode,title from ref where bibcode in ('+bibs+')', format: 'json', lang: 'adql', request :'doQuery', runid:'SimbadWeb'},
		async:false
	})
	// une fois qu'on a le resultat, remplit le titre de l'element
	.done(function(result) {
		var res = tab2map(result.data);
		if (cachebib && cachebib.length>0)
			newcache = JSON.parse(cachebib);
		$('span.bibcode').each(function(index, elemnt) {
			var bibcode = $(elemnt).text().trim();
			$(elemnt).attr('title',res[bibcode]);
			// fabrique un nouveau cache avec le nouveau titre
			newcache[bibcode]=res[bibcode];
		});
	});
	return newcache;
	

};

jQuery.fn.dataTableExt.oSort['sort-int']  = function(a,b) {
			var num;
			var str = jQuery.trim(cell.text());
			if (str.indexOf('~') >= 0 || str == '' || str == ' ' ) num = order*2147483647;
			else num = parseInt(cell.text());
			return num;
};

function disable_equinox( nb) {
	var menu = $(document.getElementsByName("frame"+nb)[0]);
	if (menu.val()=="ICRS") {
		$(document.getElementsByName("equi"+nb)[0]).attr("disabled", "disabled");
	} else {
		$(document.getElementsByName("equi"+nb)[0]).removeAttr("disabled");
	}
}

$.ajaxSetup({
	    cache: true
	 });

function showRawId() {
	if (mytable) {
		mytable.column(":contains(object count)").visible(true);
		mytable.column(":contains(object location)").visible(true);
		mytable.column(":contains(object names)").visible(true);

        $('#hideRawId').show();
        $('#showRawId').hide();
	}
};

function hideRawId() {
	if (mytable) {
		mytable.column(":contains(object count)").visible(false);
		mytable.column(":contains(object location)").visible(false);
		mytable.column(":contains(object names)").visible(false);
        $('#hideRawId').hide();

		if (mytable.column(":contains(object count)").length==0){
                $('#showRawId').hide();
        }
		else {
        	$('#showRawId').show();
		}
	}
};

// change on select[name=NbIdent]
function showSearchRawIdFixed(selected_item) {
	if (selected_item=='raw_id') {
		$('select[name=NbIdent]').parent().append('<tr id="raw_id_fixed"><td>\
			Fixed ? <input type="checkbox" id="raw_id_fixed"\
				onchange="setSearchRawIdFixed(this,true)">\
		</td></tr>');
	}
	else {
		$('input#raw_id_fixed').removeAttr('checked');
		setSearchRawIdFixed($('input#raw_id_fixed'),false);
		$('tr#raw_id_fixed').remove();
	}
}

function setSearchRawIdFixed(checkbox,selected) {
	// si coché "fixed"
	if (checkbox.checked) {
		// supprime l'ancien
		$('select[name=NbIdent]')[0].lastElementChild.remove();
		// met le nouveau
		$('select[name=NbIdent]')[0].add(new Option('raw_id','raw_id!',false,true));
	}
	else {
		// supprime l'ancien
		$('select[name=NbIdent]')[0].lastElementChild.remove();
		// met le nouveau
		$('select[name=NbIdent]')[0].add(new Option('raw_id','raw_id',false, selected));
	}
}

$.getScript("/Simlib/coo.js");
$.getScript("/Simlib/astroMath.js");
$(document).ready(function() {
	// on output options web page : clean old cookies
	if (window.location.href.indexOf("sim-fout")>0) {
	    for (var f=1 ; f <5 ; f++) {
		// Define the callback called when the 'ICRS' is selected
		var mnuOptions = $(document.getElementsByName("frame"+f)[0]);
		disable_equinox(f);
		mnuOptions.change(function(event){
		disable_equinox(event.target.name[event.target.name.length-1]);
		});
	     }
	}
	
	// convertit la table de resultat en datatable jquery
	var old = localStorage.getItem("datatable_length");
	if (old == undefined || old == "undefined") {
		old=100;
	}
	document.querySelectorAll('#datatable').forEach( function(eltDatatable) {
		mytable=        $(eltDatatable).DataTable(
			{
			 "pageLength": old, 
			 "lengthMenu": [ 100, 500, 1000, 10000 ]
			}); 
	});
	hideRawId();
	$('select[name=datatable_length]').on('change', function() {
		localStorage.setItem("datatable_length",$('select[name=datatable_length]').val());
			                });	
	//});

	// add action on change on select option raw_id for search by id page
	$('select[name=NbIdent]').on('change', function() {
		showSearchRawIdFixed($('select[name=NbIdent]').val());
			                });
	
	// min width pour que aladinlite vienne plus à gauche
	var w = $('[id=basic_data] table').width();
	$('[id=basic_data]').attr('style','min-width:'+Math.min(830,w)+'px');

	// ajoute le titre de la ref sous les bibcodes
	allbibcodes="'first'";
	cache = sessionStorage.bibcode_cache;
	if (!cache) cache="";
	$('span.bibcode').each(function(index, elemnt) {
		bib=$(elemnt).text().trim();
		if (cache.indexOf(bib)==-1) {
			if (bib.length == 19){
				allbibcodes=allbibcodes+",'"+bib+"'";
			}
		}
		else {
			var cacheBib = JSON.parse(cache);
			$(elemnt).attr('title',cacheBib[bib]);
		}
	});
	if (allbibcodes!="'first'") { 
		var newcache = getBibcodeTitles(allbibcodes,cache); 
		if (newcache)
			sessionStorage.bibcode_cache=JSON.stringify(newcache);
	}
	
	// lit le formulaire, et lance la requete sur simbad à plat, et remplit la case #nbObjects
	var showCooCount = function(coo_str) {
		$('#nbObjects').empty();
		if (!coo_str) return;
		var coo = new Coo();
		coo.prec = 10;
		coo.parse(coo_str); // + frame !
		// si on a des coordonnees et un rayon
		if ( ! isNaN(coo.lon) && ! isNaN(coo.lat) ) {
			var radiusVal = $('input[name=Radius]').val();
			var radiusUnit = $('select[name="Radius.unit"]').val();
			// Converts according units to deg
			var radius = radiusVal;
			if (radiusUnit == "arcsec") {
				radius=(radius/60)/60;
			}
			else if (radiusUnit == "arcmin") {
				radius=radius/60;
			}
			$.ajax({
				url: 'https://simbad.u-strasbg.fr/simbad/sim-tap/sync',
				type:'GET',
				cache:true,
				data: {query: 'SELECT COUNT(oid) FROM basic WHERE CONTAINS(POINT(\'ICRS\', RA, DEC), \
				CIRCLE(\'ICRS\', '+coo.lon+', '+coo.lat+', '+radius+')) = 1', format: 'json', lang: 'adql', request :'doQuery'},
				async:false
			})
			// une fois qu'on a le resultat, remplit la case du tableau avec l'info
			.done(function(result) {
				var res = result.data[0];
				$('#nbObjects').html('~ ' + res[0] + ' objects');

			});
		}

	};
	// ***** Ajoute l'info du nombre d'objets dans une recherche par coordonnées (avant de lancer sur simbad)
	// ***** grace au simbad à plat
	if ($('#sim-fcoo').length>0) {
		// ajoute dans le formulaire fcoo la case pour afficher le nombre d'objets
		$($('form[id=sim-fcoo] table tr')[5]).append('<td></td>');
		$($('form[id=sim-fcoo] table tr')[5]).append('<td id="nbObjects"></td>');
		$('#nbObjects').css('color', 'darkGreen');
		$('#nbObjects').css('font-style', 'italic');
		$('#nbObjects').attr('title','The query will return around this amount of rows');

		var showCooCountDebounced = function () {
			var coo_str = $('input[name=Coord]').val();
			debounce(showCooCount(coo_str), 400); // attend xxx ms avant de redemander
		}

		$('input[name=Coord]').on('keyup change', function() {
			showCooCountDebounced();
		});
		$('input[name=Radius]').on('keyup change', function() {
			showCooCountDebounced();
		});
		$('select[name="Radius.unit"]').on('keyup change', function() {
			showCooCountDebounced();
		});

		// initial call
		showCooCountDebounced();
	}
	if ($('#form-fid').length>0) {
		// ajoute dans le formulaire fid la case pour afficher le nombre d'objets
		$($('form[id=form-fid] table tr')[4]).append('<td></td>');
		$($('form[id=form-fid] table tr')[4]).append('<td id="nbObjects"></td>');
		$('#nbObjects').css('color', 'darkGreen');
		$('#nbObjects').css('font-style', 'italic');
		$('#nbObjects').attr('title','The query will return around this amount of rows');

		//$('#nbObjects').html('~ 0 objects');

		// lit le formulaire, et lance la requete sur simbad TAP, et remplit la case
		var showCatCount = function() {
			var name = $('input[name=Ident]').val();
			if (!name) return;

			var select = $('select[name=NbIdent]').val() ;
			// si c'est selectionne query around
			if (select == 'around') {
				$('#nbObjects').empty();
			
				// interroge Sesame
				$.ajax({
					url: 'https://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame.jsonp?',
					data: {object: name},
					method: 'GET',
					dataType: 'jsonp'
				})
				.done(function(result) {
					showCooCount(result.Target.Resolver.jpos);
				});
			}
			else {
			var cat = $('input[name=Ident]').val(); 
			if ( cat && select == 'cat' ) {
				$('#nbObjects').empty();
				$.ajax({
					url: 'https://simbad.u-strasbg.fr/simbad/sim-tap/sync',
					data: {query: 'select "size" from cat where cat_name=\''+cat+'\'', format: 'text', lang: 'adql', request :'doQuery'}
				})
				// une fois qu'on a le resultat, remplit la case du tableau avec l'info
				.done(function(result) {
					var n = result.split("\n")[2];
					if ($.isNumeric(n)) {
						$('#nbObjects').html('~ ' + n + ' objects');
					}
					else {
						console.log("erreur : reponse tap NaN : "+result);
					}

				});
			}
			}

		};

		var showCatCountDebounced = function () {
			debounce(showCatCount(), 400); // attend 200ms avant de redemander
		}

		$('input[name=Ident]').on('keyup change', function() {
			showCatCountDebounced();
		});
		$('input[name=Radius]').on('keyup change', function() {
			showCatCountDebounced();
		});
		$('select[name="NbIdent"]').on('keyup change', function() {
			showCatCountDebounced();
		});
		$('select[name="Radius.unit"]').on('keyup change', function() {
			showCatCountDebounced();
		});

		// initial call
		showCatCountDebounced();
	}
	else if ($('form[action=sim-ref]').length>0 ) { // id=sim-fref
               // ajoute dans le formulaire la case pour afficher le nombre d'objets
                $($('tr:has([name=simbo]) td')[0]).attr('id','nbObjects');
                $('#nbObjects').css('color', 'darkGreen');
                $('#nbObjects').css('font-style', 'italic');
                $('#nbObjects').css('text-align','right');
		$('#nbObjects').attr('title','The query will return around this amount of rows');

                //$('#nbObjects').html('~ 0 objects');

                // lit le formulaire, et lance la requete sur simbad TAP, et remplit la case
                var showBibCount = function() {
                        var bib = $('input[name=bibcode]').val();
			if (bib.startsWith("10.")) {field="doi";} else {field="bibcode";}
			if (!bib) {
				$('#nbObjects').empty();
				return;
			}
		
			$.ajax({
				url: 'https://simbad.u-strasbg.fr/simbad/sim-tap/sync',
				data: {query: 'select nbobject from ref where '+field+'=\''+bib+'\'', 
					format: 'text', lang: 'adql', request :'doQuery'}
			})
			// une fois qu'on a le resultat, remplit la case du tableau avec l'info
			.done(function(result) {
				$('#nbObjects').empty();
				var n = result.split("\n")[2];
				if ($.isNumeric(n)) {
					$('#nbObjects').html('~ ' + n + ' objects');
				}
				else {
					console.log("erreur : reponse tap NaN : "+result);
				}

			});
		};
               var showBibCountDebounced = function () {
                        debounce(showBibCount(), 400); // attend 200ms avant de redemander
                }

                $('input[name=bibcode]').on('keyup change', function() {
                        showBibCountDebounced();
                });
		// initial call
		showBibCountDebounced();

	
	}

});


/*
$(document).ready(function() {
$('[name=Ident]').autocomplete({
      delay: 300,
      minLength: 2,
      source: function(request, response) {
          $.getJSON("http://simbad.u-strasbg.fr/tools/suggestNames/search", {
              kw: request.term,
          }, function(result) {
              var array = result.status!="success" ? [] : $.map(result.data, function(item) {
                  return {
                      label: item.objName + ' (' + item.nbRef + ' refs)',
                      value: item.objName,
                      oid: item.oid
                  };
              });
              response(array);
          });
      },
      focus: function(event, ui) {
                    // prevent autocomplete from updating the textbox
                    event.preventDefault();
                }//,
     // To add an additional action when user selects a suggestion (default action: copy value in text field)
     //  select: function(event, ui) {
     //           }
    });
  });
  */
