var aladin=null;
var mySrc = [];
var catObjects = null;
var divTypes = null;
var keyOpen = "wantedAL"+document.URL.substring(document.URL.lastIndexOf("/"),document.URL.indexOf("?"));
var showRadius=false;
var useStorage=false;
var searchCircle = null;
var simbadCatalog=null;
var simbadCatalogStick=false;
var targetSimbad=null;

/* fonction pour eviter de lancer trop souvent une méthode
 * utilisée quand on veut faire un appel à chaque touche */
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

var changeALDebounced = function () {
	var zoom ="";
	var coo="";
	  if ($('input[name=Coord]')) {
		  coo = $('input[name=Coord]').val();
	  }
	zoom=readRadius();
	debounce(moveAladinLitePreview(coo,zoom), 600); // attend xxx ms avant de redemander
};

/**
* Fonction qui extrait du formulaire des coordonnées le rayon a la bonne unité
*/
var readRadius = function () {
	var zoom = null;
	  if ($('input[name=Radius]')[0]) {
		  unit= $('[name="Radius.unit"]').val();
		  value= $('[name="Radius"]').val();
		  if (value) {
			  if (unit=="arcmin") {
				  zoom = $('input[name=Radius]').val()/60;
			  } else if (unit=="deg") {
				  zoom = $('input[name=Radius]').val();
			  } else { // sec
				zoom = $('input[name=Radius]').val()/(60*60);
			  }
		  }
	  }
	if (zoom==null) zoom=0.01;
	return zoom;
}


$(document).ready(function() {
	if(typeof(Storage) !== "undefined") { useStorage=true; }
	$('input[name=Coord]').on('keyup change', function() {
		if (aladin) { changeALDebounced();
		createSimbadCenter(coo.lon, coo.lat); }
	});
	$('[name^=Radius]').on('keyup change', function() {
		if (aladin) changeALDebounced();
	});
	// form "Query by Coordinates" -> always show the circle
	if ($('#sim-fcoo')[0] != undefined) {
		showRadius=true;
	}
	else { // for other forms, just show the circle if the mouse is in the radius form
		if ($('input[name=Radius]')[0]) {
			$('input[name=Radius]').mouseenter(function() {
			  newRadius=$('input[name=Radius]').val();
			  if (newRadius!="" && $.isNumeric(newRadius)) {
				showRadius=true;
				searchCircle.setRadius((newRadius/60));
				searchCircle.show();
			  }
			  else {
				searchCircle.hide();
			  }
			});
			$('input[name=Radius]').mouseleave(function() {
			  searchCircle.hide();
			});

		}
	}
    
	const eltDatatable = document.querySelector('#datatable.astrobject-list');
	if (eltDatatable){
		// Index all lines for safe object search:
		eltDatatable.querySelectorAll('tr')
		.forEach(line => {
			const firstCol = line.querySelector('td:first-child');
			if (firstCol)
				line.setAttribute('data-idx', firstCol.innerText.trim());
		});

		/* Split the table area between the data table and AladinLite
		 * (each part being resizable): */
		initSplitTableAladin(eltDatatable);

		// Initialize AladinLite:
		const coordstring = eltDatatable.getAttribute('data-coordstring');
		const radius      = Number.parseFloat(eltDatatable.getAttribute('data-radius'));
		const survey      = eltDatatable.getAttribute('data-survey');
		myExtendedAladinLiteList(survey, (coordstring ? coordstring : ""), radius);
	}

	 const plot_button = document.querySelector('input[name="plot"]');
        if (plot_button) {
                plot_button.setAttribute("onclick"," \
					parentElement.parentElement.querySelector('.divided .divider').dispatchEvent (\
						new MouseEvent('mousedown', { clientX: document.documentElement.clientWidth})\
					);\
					window.dispatchEvent (\
						new MouseEvent('mousemove', { clientX: document.documentElement.clientWidth*0.7})\
				  	);\
					window.dispatchEvent (new MouseEvent('mouseup'));\
			");

		}
});

//  $("head").append('<link rel="stylesheet" href="http://ajax.googleapis.com/ajax/libs/jqueryui/1.11.24/themes/smoothness/jquery-ui.css" />');
$("head").append('<link rel="stylesheet" href="//simbad.u-strasbg.fr/Simlib/jquery/jquery-ui.css" />');

var cachedScript = 
function( url, options ) {
	// Allow user to set any option except for dataType, cache, and url
	options = $.extend( options || {}, {
		dataType: "script",
		cache: true,
		url: url
	});
	// Use $.ajax() since it is more flexible than $.getScript
	// Return the jqXHR object so we can chain callbacks
	return $.ajax( options );
};

function myAladinLite(s, id, zoom, setsize=true) {
$("head").append("<link rel='stylesheet' type='text/css' href='//aladin.u-strasbg.fr/AladinLite/api/v2/latest/aladin.min.css' />");
//$.getScript("http://ajax.googleapis.com/ajax/libs/jqueryui/1.8.24/jquery-ui.min.js", function() {
        cachedScript("//simbad.u-strasbg.fr/Simlib/aladin.min.js").done(function() {
	// set size of aladin-lite div
	if (setsize)	$('#aladin-lite-div').css({'width':'260px','height':'250px'});

	// extract Simbad coordinates from the web page
	var coo=null;
	var center = $('#Coord');
    if (center.length==0) center = $('input[name=Coord]');

	if (center.length > 0) {
		coo = new Coo();
		id = center.val();
		if (id == "") id ="0+0"; 
		if (!coo.parse(id)) {
			id ="0+0";coo.parse(id); }
//		window.console.log(center.val());
    		initSurveys(coo.lon, coo.lat);
	}

	if (zoom<0.01) {
		zoom=0.01;
		console.log("zoom too small to display images => reset to 0.01d");
	}
        aladin = A.aladin('#aladin-lite-div', {survey:s, zoom:zoom, target:id , showZoomControl:false, showGotoControl:false, showFrame:false, showLayersControl:false, showReticle:false, showSimbadPointerControl:true }); //reticleColor:"#30a090", reticleSize:30});
	$('div.aladin-simbadPointerControl-container').css({'top':'200px'});
	$('.aladin-simbadPointerControl').css({'width':'22px','height':'22px'});
	$('div.aladin-logo-container a').attr("href",'https://aladin.cds.unistra.fr/AladinLite/?target='+encodeURIComponent(id)+'&fov='+zoom*2+'&survey='+encodeURIComponent(s));
        selectSurvey(s);
        aladin.setFOVRange(0.01,103);

	addSimbadCone();
	// create a target overlayed for the Simbad coordinates 
	createSimbadCenter(coo.lon, coo.lat);
	if (simbadCatalogStick) showSimbadCone();

        // prepare overlays
        var overlay = A.graphicOverlay({color: 'cyan'});
        aladin.addOverlay(overlay);
	
        searchCircle = A.circle(coo.lon, coo.lat, readRadius());
        searchCircle.hide();
        overlay.add(searchCircle);
	if (showRadius) searchCircle.show();
        
        // radius given in the SED form 
	if ($('#FormSED')[0] != undefined ) {
		// add listener when mouse enters/leaves query around form
		$('#FormSED').mouseenter(function() {
			newRadius=$('#radiusSED').val();
		  if (newRadius!="" && $.isNumeric(newRadius)) {
			showRadius=true;
		      searchCircle.setRadius((newRadius/60)/60);
			  searchCircle.show();
			}
		});
		$('#FormSED').mouseleave(function() {
		  searchCircle.hide();
		});
		
		// add listener when value of radius is updated
		$('#radiusSED').bind("change keyup", function() {
		  newRadius = $(this).val();
		  if (newRadius=="" || ! $.isNumeric(newRadius)) {
		      searchCircle.hide();
			showRadius=false;
		      return;
		    }
		  else {
		      searchCircle.setRadius((newRadius/60)/60);
			showRadius=true;
		      searchCircle.show();
		  }
		});
	}       
	if (showRadius) searchCircle.show();
});
}

function myAladinLiteList(s, center, zoom) {
	$("head").append("<link rel='stylesheet' type='text/css' href='//aladin.u-strasbg.fr/AladinLite/api/v2/latest/aladin.min.css' />");
	//$("head").append('<link rel="stylesheet" href="http://ajax.googleapis.com/ajax/libs/jqueryui/1.8.24/themes/smoothness/jquery-ui.css" />');

	cachedScript("//simbad.u-strasbg.fr/Simlib/jquery-ui.min.js").done(function() {
		//cachedScript("http://ajax.googleapis.com/ajax/libs/jqueryui/1.8.24/jquery-ui.min.js").done(function() {
		cachedScript("//aladin.u-strasbg.fr/AladinLite/api/v2/latest/aladin.js").done(function() {

		//cachedScript("//simbad.u-strasbg.fr/Simlib/aladin.min.js").done(function() {
			myjquery = $;
			aladin = A.aladin('#aladin-lite-div', {survey:s, zoom:zoom, target:center , showZoomControl:false, showGotoControl:false, showFrame:false, showLayersControl:false, reticleColor:"#30a090", reticleSize:30, showReticle:false });

			selectSurvey(s);
			aladin.setFOVRange(0.003,103);
			aladinLiteCatalog();
			linkSimbadWithAladinLite();
			var coo = new Coo();
			coo.parse(center);
			createSimbadCenter(coo.lon, coo.lat);

			// on selectionne la première ligne
			var premiereLigne = $('#datatable tr').eq(1);
			selectThisLineInAladin(premiereLigne);
			selectLine(premiereLigne);
			attachCloseEvent();
			/*
				if(isAladinLiteWanted()) {
					openAladinLiteDialog();
				}
				else {
					initAladinLiteDialog();
					$("#aladinLiteDialog").dialog("close");
				}*/
		});
	});
}

function myExtendedAladinLiteList(s, center, zoom){
	// Load styles of AladinLite:
	$("head").append("<link rel='stylesheet' type='text/css' href='//aladin.u-strasbg.fr/AladinLite/api/v2/latest/aladin.min.css' />");
	
	// check empty information -> set to default
	if (! center || center.trim().length==0) {
		center=getRowCooStr($('#datatable tr').eq(1));
	}
	if (! zoom || isNaN(zoom)) {
		zoom=readRadius();
	}
	if (! s || s.trim().length==0) {
		s=Object.keys(surveyId2Name)[1];
	}

	// Load AladinLite:
	cachedScript("//aladin.u-strasbg.fr/AladinLite/api/v2/latest/aladin.js")
	.done(() => initExtendedAladinLiteList(s, center, zoom));
}

function initExtendedAladinLiteList(s, center, zoom){
	let drawPM  = true;
	let pmYears = 1000; // show PM after pmYears years

	const SUPPORTED_TYPES = {
		star  : { label: "Stars"           , example: "*"  , test: otype=>/^(\*|\*iC|\*iN|\*iA|V\*\?|Pe\*|HB\*|Y\*O|Ae\*|Em\*|Be\*|BS\*|RG\*|AB\*|C\*|S\*|sg\*|s\*r|s\*y|s\*b|HS\*|pA\*|WD\*|LM\*|BD\*|N\*|OH\*|pr\*|TT\*|WR\*|PM\*|HV\*|V\*|Ir\*|Or\*|Er\*|RC\*|RC\?|Ro\*|a2\*|Psr|BY\*|RS\*|Pu\*|RR\*|Ce\*|dS\*|RV\*|WV\*|bC\*|cC\*|gD\*|SX\*|LP\*|Mi\*|SN\*|su\*|Pl\?|Pl)$/.test(otype)},
		galaxy: { label: "Galaxies"        , example: "G"  , test: otype=>/^(G|PoG|GiC|BiC|GiG|GiP|HzG|rG|H2G|LSB|AG\?|Q\?|Bz\?|BL\?|EmG|SBG|bCG|LeI|LeG|LeQ|AGN|LIN|SyG|Sy1|Sy2|Bla|BLL|OVV|QSO)$/.test(otype)},
		nebula: { label: "Nebula, PN, SNR" , example: "PN" , test: otype=>/^(PN\?|PN|GNe|BNe|DNe|RNe|SR\?|SNR)$/.test(otype)},
		hii   : { label: "HII regions"     , example: "HII", test: otype=>/^(HII|sh)$/.test(otype)},
		radio : { label: "Radio, HI, Maser", example: "Rad", test: otype=>/^(Rad|mR|cm|mm|smm|HI|rB|Mas)$/.test(otype)},
		ir    : { label: "IR objects"      , example: "IR" , test: otype=>/^(IR|FIR|MIR|NIR)$/.test(otype)},
		uv    : { label: "UV objects"      , example: "UV" , test: otype=>"UV"===otype},
		xray  : { label: "X-ray objects"   , example: "X"  , test: otype=>/^(X|UX\?|ULX)$/.test(otype)}
	};

	// Create the AladinLite view:
	aladin = A.aladin('#aladin-lite-div', {
		survey           : s,
		zoom             : zoom,
		target           : center,
		showZoomControl  : false,
		showGotoControl  : false,
		showFrame        : false,
		showShareControl : false
	});

	// Select the desired image survey:
	selectSurvey(s);

	// Adapt the Field Of View:
	aladin.setFOVRange(0.003,103);

	// Create a catalogue with all objects of the datatable:
	catObjects = aladinLiteCatalog(drawFunction, false);

	// Set the legend into AladinLite:
	const formerLegend = document.getElementById('aladin-lite-div').querySelector('.aladin-legend-container');
	if (formerLegend)
		formerLegend.remove();
	addLegend();

	// Add sources and update the legend in AladinLite:
	updateAladinLiteCatalog();

	// Link each datatable's row to the corresponding object in AladinLite:
	linkSimbadWithAladinLite();

	// Center the view on the queried object (i.e. given coordinates):
	var coo = new Coo();
	coo.parse(center);
	createSimbadCenter(coo.lon, coo.lat);

	/* Update the catalog and the legend whenever the object set changes:
	 * (ex: change page, search, sort column, ...) */
	$('#datatable').on('draw.dt', updateAladinLiteCatalog);

	// By default, select the first row:
	var premiereLigne = $('#datatable tr').eq(1);
	selectThisLineInAladin(premiereLigne);
	selectLine(premiereLigne);

	function updateAladinLiteCatalog(){
		// Update the sources visible in AladinLite:
		refreshAladinLiteCatalog(catObjects);

		// Update the AladinLite's legend:
		refreshLegend();
	}

	function parseOtype(otype){
		let category = "other";
		Object.keys(SUPPORTED_TYPES).forEach(t=>{
			if (SUPPORTED_TYPES[t].test(otype)){
				category = t;
				return;
			}
		});
		return category;
	}

	function drawFunction(source, canvasCtx) {

		sourceSize = 6;
		
		canvasCtx.beginPath();
	
		canvasCtx.strokeStyle = 'red';
		canvasCtx.lineWidth   = 2;
		canvasCtx.globalAlpha = 1;
		canvasCtx.lineJoin    = 'round';
		canvasCtx.lineCap     = 'round';
		canvasCtx.fillStyle   = '#FFF1';
	
		// Resolve the object type:
		const optType = parseOtype(source.data.otype);

		// If no visibility option for this object type, consider it as TRUE:
		if (!aladin.options.simbad_types[optType])
			aladin.options.simbad_types[optType] = { show: true };
		
		// Adapt the symbol in function of the object type:
		switch(optType){
			// CASE: Star => red cercle
			case 'star':
				if (aladin.options.simbad_types[optType].show){
					canvasCtx.arc(source.x, source.y, sourceSize, 0, 2 * Math.PI, false);
					canvasCtx.stroke();
					canvasCtx.fill();
				}
				break;
	
			// CASE: Galaxy => blue ellipse
			case 'galaxy':
				if (aladin.options.simbad_types[optType].show){
					canvasCtx.strokeStyle = 'white';
					canvasCtx.lineWidth   = 4;
					canvasCtx.globalAlpha = .5;
					canvasCtx.ellipse(source.x, source.y, sourceSize*1.5, sourceSize, (source.data.galdim_angle ? (3*Math.PI/2.)-source.data.galdim_angle*Math.PI/180. : 0), 0, 2 * Math.PI);
					canvasCtx.stroke();
					canvasCtx.fill();
					canvasCtx.lineWidth   = 2;
					canvasCtx.globalAlpha = 1;
					canvasCtx.strokeStyle = 'blue';
					canvasCtx.ellipse(source.x, source.y, sourceSize*1.5, sourceSize, (source.data.galdim_angle ? (3*Math.PI/2.)-source.data.galdim_angle*Math.PI/180. : 0), 0, 2 * Math.PI);
					canvasCtx.stroke();
					canvasCtx.fill();
				}
				break;
	
			// CASE: Nebula, PN, SNR => red square
			case 'nebula':
				if (aladin.options.simbad_types[optType].show){
					canvasCtx.moveTo(source.x-sourceSize, source.y-sourceSize);
					canvasCtx.lineTo(source.x+sourceSize,  source.y-sourceSize);
					canvasCtx.lineTo(source.x+sourceSize,  source.y+sourceSize);
					canvasCtx.lineTo(source.x-sourceSize,  source.y+sourceSize);
					canvasCtx.lineTo(source.x-sourceSize, source.y-sourceSize);
					canvasCtx.stroke();
					canvasCtx.fill();
				}
				break;
	
			// CASE: HII region, HIshell => small red circle
			case 'hii':
				if (aladin.options.simbad_types[optType].show){
					canvasCtx.moveTo(source.x-sourceSize, source.y);
					canvasCtx.lineTo(source.x+sourceSize, source.y);
					canvasCtx.stroke();
					
					canvasCtx.moveTo(source.x, source.y-sourceSize);
					canvasCtx.lineTo(source.x, source.y+sourceSize);
					canvasCtx.stroke();
				}
				break;
	
			// CASE: Radio, HI, Maser => blue triangle
			case 'radio':
				if (aladin.options.simbad_types[optType].show){
					canvasCtx.strokeStyle = 'white';
					canvasCtx.lineWidth   = 4;
					canvasCtx.globalAlpha = .5;
					canvasCtx.moveTo(source.x, source.y-sourceSize);
					canvasCtx.lineTo(source.x+sourceSize, source.y+sourceSize);
					canvasCtx.lineTo(source.x-sourceSize, source.y+sourceSize);
					canvasCtx.lineTo(source.x, source.y-sourceSize);
					canvasCtx.stroke();
					canvasCtx.fill();
					canvasCtx.strokeStyle = 'blue';
					canvasCtx.lineWidth   = 2;
					canvasCtx.globalAlpha = 1;
					canvasCtx.moveTo(source.x, source.y-sourceSize);
					canvasCtx.lineTo(source.x+sourceSize, source.y+sourceSize);
					canvasCtx.lineTo(source.x-sourceSize, source.y+sourceSize);
					canvasCtx.lineTo(source.x, source.y-sourceSize);
					canvasCtx.stroke();
					canvasCtx.fill();
				}
				break;
	
			// CASE: IR => red losange
			case 'ir':
				if (aladin.options.simbad_types[optType].show){
					canvasCtx.moveTo(source.x, source.y-sourceSize);
					canvasCtx.lineTo(source.x+sourceSize, source.y);
					canvasCtx.lineTo(source.x, source.y+sourceSize);
					canvasCtx.lineTo(source.x-sourceSize, source.y);
					canvasCtx.lineTo(source.x, source.y-sourceSize);
					canvasCtx.stroke();
					canvasCtx.fill();
				}
				break;
	
			// CASE: UV => pink star
			case 'uv':
				if (aladin.options.simbad_types[optType].show){
					canvasCtx.strokeStyle = '#f800d2';
					canvasCtx.lineWidth = 1;
					canvasCtx.moveTo(source.x, source.y-sourceSize);
					canvasCtx.lineTo(source.x, source.y+sourceSize);
					canvasCtx.stroke();
					canvasCtx.moveTo(source.x-sourceSize, source.y);
					canvasCtx.lineTo(source.x+sourceSize, source.y);
					canvasCtx.stroke();
					sourceSize = sourceSize*0.8;
					canvasCtx.moveTo(source.x-sourceSize, source.y-sourceSize);
					canvasCtx.lineTo(source.x+sourceSize, source.y+sourceSize);
					canvasCtx.stroke();
					canvasCtx.moveTo(source.x+sourceSize, source.y-sourceSize);
					canvasCtx.lineTo(source.x-sourceSize, source.y+sourceSize);
					canvasCtx.stroke();
				}
				break;
	
			// CASE: X-Ray => black cross
			case 'xray':
				if (aladin.options.simbad_types[optType].show){
					sourceSize = sourceSize*0.8;
					canvasCtx.strokeStyle = 'yellow';
					canvasCtx.lineWidth   = 4;
					canvasCtx.globalAlpha = .5;
					canvasCtx.moveTo(source.x-sourceSize, source.y-sourceSize);
					canvasCtx.lineTo(source.x+sourceSize, source.y+sourceSize);
					canvasCtx.stroke();
					canvasCtx.moveTo(source.x+sourceSize, source.y-sourceSize);
					canvasCtx.lineTo(source.x-sourceSize, source.y+sourceSize);
					canvasCtx.stroke();
					canvasCtx.strokeStyle = 'black';
					canvasCtx.lineWidth   = 2;
					canvasCtx.globalAlpha = 1;
					canvasCtx.moveTo(source.x-sourceSize, source.y-sourceSize);
					canvasCtx.lineTo(source.x+sourceSize, source.y+sourceSize);
					canvasCtx.stroke();
					canvasCtx.moveTo(source.x+sourceSize, source.y-sourceSize);
					canvasCtx.lineTo(source.x-sourceSize, source.y+sourceSize);
					canvasCtx.stroke();
				}
				break;
	
			// DEFAULT: => yellow circle
			default:
				if (aladin.options.simbad_types[optType].show){
					canvasCtx.strokeStyle = 'black';
					canvasCtx.lineWidth   = 4;
					canvasCtx.globalAlpha = .5;
					canvasCtx.arc(source.x, source.y, sourceSize, 0, 2 * Math.PI, false);
					canvasCtx.stroke();
					canvasCtx.strokeStyle = 'yellow';
					canvasCtx.lineWidth   = 2;
					canvasCtx.globalAlpha = 1;
					canvasCtx.arc(source.x, source.y, sourceSize, 0, 2 * Math.PI, false);
					canvasCtx.stroke();
				}
				break;
		}

		// add proper motion symbol (i.e. an arrow):
		addPM(source, canvasCtx);
		
		canvasCtx.closePath();
	}

	function addLegend(){
		const divAladin = document.getElementById('aladin-lite-div');
	
		// Create the legend container:
		const divContainer = document.createElement('div');
		addClass(divContainer, 'aladin-legend-container');
		divAladin.appendChild(divContainer);
	
		// Create the legend box:
		const divLegend = document.createElement('div')
		addClass(divLegend, 'aladin-legend');
		addClass(divLegend, 'aladin-box');
		divContainer.appendChild(divLegend);
		// Close button:
		const btnClose = document.createElement('a');
		addClass(btnClose, 'aladin-closeBtn');
		btnClose.innerHTML = 'x';
		addEvent(btnClose, 'click', ()=>removeClass(divContainer, 'is-active'));
		divLegend.appendChild(btnClose);
		// Title:
		const divTitle = document.createElement('div');
		addClass(divTitle, 'aladin-label');
		divTitle.innerHTML = "Legend";
		divLegend.appendChild(divTitle);
	
		// Append a line for all symbols:
		divTypes = document.createElement('div');
		addClass(divTypes, 'symbol-list');
		divLegend.appendChild(divTypes);
	
		// Append a line for options:
		divLegend.appendChild(document.createElement('hr'));
		// ...draw proper motion:
		const divPMLine = document.createElement('div');
		addClass(divPMLine, 'symbol-line');
		const chkPM = document.createElement('input');
		chkPM.setAttribute('type', 'checkbox');
		chkPM.setAttribute('id', 'chkPM');
		if (drawPM)
			chkPM.setAttribute('checked', 'checked');
		chkPM.addEventListener('change', ()=>{
			drawPM = !drawPM;
			aladin.view.requestRedraw();
		});
		divPMLine.appendChild(chkPM);
		const lblPM = document.createElement('label');
		lblPM.setAttribute('for', 'chkPM');
		lblPM.innerText = 'Show proper motions';
		divPMLine.appendChild(lblPM);
		divLegend.appendChild(divPMLine);
		// ...proper motion time:
		const divPMTime = document.createElement('div');
		addClass(divPMTime, 'symbol-line');
		const rgPM = document.createElement('input');
		rgPM.setAttribute('id', 'pmTime');
		rgPM.setAttribute('type', 'range');
		rgPM.setAttribute('min', '100');
		rgPM.setAttribute('max', '100000');
		rgPM.setAttribute('step', '100');
		rgPM.setAttribute('value', pmYears);
		divPMTime.appendChild(rgPM);
		const lblPMTime = document.createElement('label')
		lblPMTime.setAttribute('for', 'pmTime');
		lblPMTime.innerText = 'PM time';
		divPMTime.appendChild(lblPMTime);
		const spanTime = document.createElement('span');
		spanTime.innerText = ` (${pmYears} yr)`;
		lblPMTime.appendChild(spanTime);
		rgPM.addEventListener('change', ()=>{
			spanTime.innerText = ` (${rgPM.value} yr)`;
			pmYears = rgPM.value;
			aladin.view.requestRedraw();
		});
		rgPM.addEventListener('input', ()=>{
			spanTime.innerText = ` (${rgPM.value} yr)`;
			pmYears = rgPM.value;
			aladin.view.requestRedraw();
		});
		divLegend.appendChild(divPMTime);
	
		// Finally add the trigger button (showing the legend box):
		const btnLegend = document.createElement('div');
		addClass(btnLegend, 'aladin-legend-button');
		btnLegend.setAttribute('title', "Show legend");
		addEvent(btnLegend, 'click', ()=>addClass(divContainer, 'is-active'));
		divContainer.appendChild(btnLegend);
	}

	function refreshLegend(){
		// Nothing to do, if no legend yet:
		if (divTypes){
			// Remove previous legend item:
			divTypes.innerHTML = "";

			// Get all extracted object types:
			const otypes = catObjects.sources.map(s => s.data.otype);

			// Count objects by type:
			aladin.options.simbad_types = otypes.map(obj=>parseOtype(obj))
											.reduce((accu, type)=>{
												if (accu[type])
													accu[type].count++;
												else
													accu[type] = { show: true, count: 1 };
												return accu;
											}, {});

			// Insert a new line in the legend box for each found object type:
			Object.entries(aladin.options.simbad_types)
				.sort((t1,t2)=>t2[1].count-t1[1].count)
				.forEach(t=>{
				if (t[0] === "other")
					divTypes.appendChild(buildSymbolLine(t[0], t[0], "Other types", t[1].count));
				else{
					const typeDesc = SUPPORTED_TYPES[t[0]];
					divTypes.appendChild(buildSymbolLine(typeDesc.example, t[0], typeDesc.label, t[1].count));
				}
				}
			);
		}
	}
	
	function buildSymbolLine(type, optName, text, count){
		const divLine = document.createElement('div');
		addClass(divLine, 'symbol-line');
		
		const input = document.createElement('input');
		input.setAttribute('type', "checkbox");
		const id = `symbol_${type}`;
		input.setAttribute('id', id);
		input.setAttribute('checked', "checked");
		addEvent(input, 'change', ()=>{
			aladin.options.simbad_types[optName].show = !aladin.options.simbad_types[optName].show;
			catObjects.sources.forEach(s => {
				s.isShowing = aladin.options.simbad_types[parseOtype(s.data.otype)].show;
			});
			aladin.view.requestRedraw();
		});
		divLine.appendChild(input);
	
		const label = document.createElement('label');
		label.setAttribute('for', id);
		divLine.appendChild(label);
	
		const canvas = document.createElement('canvas');
		canvas.setAttribute('width' , 20);
		canvas.setAttribute('height', 20);
		drawFunction({x: 10, y: 10, data: {"otype": type}}, canvas.getContext("2d"));
		label.appendChild(canvas);
	
		label.appendChild(document.createTextNode(`${text} (${(count ? count : '?')})`))
	
		return divLine;
	}

	function addPM(source, canvasCtx){

		if (drawPM && source.data.pmra && source.data.pmdec){
			
			/* Compute the point corresponding to the extremety of the arrow
			 * representing the proper motion: */
			const pm_coords = applyPM(source, aladin.view, sourceSize);
	
			if (pm_coords){
			
				// Compute vectors, angle and length:
				// ...vector for the ordinate axis (relative to the target source):
				const S = [ 1 , 0 ];
				// ...vector for the proper motion (relative to the target source):
				const P = [ pm_coords.x-source.x , pm_coords.y-source.y ]; //
				// ...angle between both vectors:
				const angle  = Math.acos(((S[0]*P[0])+(S[1]*P[1]))/(Math.sqrt(S[0]*S[0]+S[1]*S[1])*Math.sqrt(P[0]*P[0]+P[1]*P[1])));
				// ...length of the proper motion vector:
				const length = Math.sqrt( P[0]*P[0] + P[1]*P[1] );
	
				const arrowWidth = 1;
	
				if (length > sourceSize+arrowWidth){
	
					// Draw the proper motion arrow's line:
					canvasCtx.lineWidth = 1;
					canvasCtx.moveTo(source.x, source.y);
					canvasCtx.lineTo(pm_coords.x, pm_coords.y);
					canvasCtx.stroke();
	
					// Compute coordinates of the arrow elements (for the moment on the ordinate axis):
					let arrow_up  = [ S[0]+length-sourceSize , S[1]-sourceSize/2. ];
					let arrow_dwn = [ S[0]+length-sourceSize , S[1]+sourceSize/2. ];
	
					// Rotate these arrow elements to align with the arrow's line:
					arrow_up  = [ arrow_up[0]*Math.cos(angle)-arrow_up[1]*Math.sin(angle)   , arrow_up[0]*Math.sin(angle)+arrow_up[1]*Math.cos(angle)];
					arrow_dwn = [ arrow_dwn[0]*Math.cos(angle)-arrow_dwn[1]*Math.sin(angle) , arrow_dwn[0]*Math.sin(angle)+arrow_dwn[1]*Math.cos(angle)];
	
					// Draw the arrow elements:
					canvasCtx.moveTo(pm_coords.x, pm_coords.y);
					canvasCtx.lineTo(source.x+arrow_up[0], source.y+(P[1] < 0 ? -1 : 1)*arrow_up[1]);
					canvasCtx.stroke();
					canvasCtx.moveTo(pm_coords.x, pm_coords.y);
					canvasCtx.lineTo(source.x+arrow_dwn[0], source.y+(P[1] < 0 ? -1 : 1)*arrow_dwn[1]);
					canvasCtx.stroke();
				}
			}
		}
	}
	
	function applyPM(source, view, sourceSize){
	
		// Convert Proper Motion (mas/yr) into degrees:
		const coords = {
			ra: source.ra+source.data.pmra*pmYears/3600000.,
			dec: source.dec+source.data.pmdec*pmYears/3600000.
		};
	
		let xy;
		if (view.cooFrame.system != CooFrameEnum.SYSTEMS.J2000) {
			const lonlat = CooConversion.J2000ToGalactic([coords.ra, coords.dec]);
			xy = view.projection.project(lonlat[0], lonlat[1]);
		}
		else {
			xy = view.projection.project(coords.ra, coords.dec);
		}
	
		if (xy) {
			const xyview = AladinUtils.xyToView(xy.X, xy.Y, view.width, view.height, view.largestDim, view.zoomFactor, true);
			if (xyview) {
				// check if source is visible in view
				if (xyview.vx<=(view.width+sourceSize) && xyview.vx>=(0-sourceSize) &&
					xyview.vy<=(view.height+sourceSize) && xyview.vy>=(0-sourceSize)
					&& (source.x != xyview.vx || source.y != xyview.vy)) {
					return { x: xyview.vx, y: xyview.vy };
				}
			}
		}
	
		return null;
	}
}

/**
 * Fonction qui cree un catalogue pour le point central
 */ 
function createSimbadCenter(ra, dec) {
	if (targetSimbad) {
		targetSimbad.removeAll();
                aladin.view.requestRedraw();
                targetSimbad=null;
	}

	targetSimbad = A.catalog({color: '#30a090', shape: 'plus', sourceSize: 18});
	aladin.addCatalog(targetSimbad);
	targetSimbad.addSources([A.source(ra, dec)]);
}

/**
 * Fonction qui va lier la fermeture de la dialog JQueryUI avec une portion de code
 */
function attachCloseEvent() {
  $("#aladinLiteDialog").bind("dialogclose", function(event) {
    $(".imgpreview").css("visibility", "hidden");
    $("#datatable tr").removeClass("selectedRow");
    $("#datatable tr").removeClass("hoverRow");
    if (useStorage) localStorage.setItem(keyOpen, "false");
  });
}

/**
 * Fonction qui dit si oui ou non l'utilisateur veut de la fenetre AladinLite
 * @return {Boolean} Vrai si l'utilisateur en veut, faux sinon
 */
function isAladinLiteWanted() {
  var res = true;
  if ($('input[name=coo]').length == 0) {// si on n'a pas un champ coordonnees
	  res = false;
  }
  if(useStorage)  {
	var old = localStorage.getItem(keyOpen);
	if (old) {
		res = (old === "true");
	}
	// else si on a n'a pas deja enregistre l'info on reste avec le defaut
  }

  return res;
}

/**
 * Extrait une coordonnée (string) depuis une ligne du tableau
 */
function getRowCooStr(ligne) {
  var index_col_ra = getIndexOf("RA");
  var index_col_dec = getIndexOf("DEC");

  var raStr = ligne.find('td').eq(index_col_ra).html().trim();
  var decStr = ligne.find('td').eq(index_col_dec).html().trim();

  return raStr + " " + decStr;
}

/**
 * Fonction qui prend en entrée une ligne et qui positionne AladinLite
 * sur l'objet correspondant à cette ligne de Simbad
 * @param  {JQueryObject} ligne La ligne Simbad sur laquelle positionner AladinLite
 */
function selectThisLineInAladin(ligne) {
  var coo = new Coo();
  coo.parse(getRowCooStr(ligne));
  aladin.gotoPosition(coo.lon, coo.lat);
  return coo;
}

/**
 * Fonction qui permet de selectionner la ligne passée en paramêtre
 * @param  {ObjectJQuery} line La ligne a selectionner
 */
function selectLine(line) {
  var images = $(".imgpreview");
  var image_sel = line.find("img[class='imgpreview']");

  images.css("visibility", "hidden");
//  image_sel.css("visibility", "visible");

  $("#datatable tr").removeClass("selectedRow");
  line.addClass("selectedRow");
}

/**
 * Fonction qui permet de créer un catalogue dans aladinLite à partir des
 * enregistrements présent dans la page Simbad
 * 
 * @param {Function} [drawFunction] Fonction a utiliser pour dessiner chaque
 *                                  source du catalogue a creer.
 *                                  Si omit, forme par defaut.
 * @param {boolean} [withInit] Si TRUE, la fonction initialisera le
 *                             catalogue avec les premieres lignes du DataTable,
 *                             mais si FALSE, le catalogue sera cree vide.
 * 
 * @returns {Object} Le nouveau catalogue.
 */
function aladinLiteCatalog(drawFunction, withInit=true)
{
	// define function triggered when an object is hovered
	aladin.on('objectHovered', function(object) {
		//catalog.deselectAll();
		//object.select();
		if (object) {
			$('#datatable tr').removeClass('hoverRow');
			$(`#datatable tr[data-idx="${object.data.idx}"]`).addClass('hoverRow');
		}
		else {
			$('#datatable tr').removeClass('hoverRow');
		}
	});

	// define function triggered when an object is clicked
	aladin.on('objectClicked', function(object) {
		if (object['data']) {
			const ligne = $(`#datatable tr[data-idx="${object.data.idx}"]`);

			// on met aladinLite à la position du click :
			selectThisLineInAladin(ligne);

			// on selectionne la ligne sur l'interface Simbad :
			selectLine(ligne);

			// scrolling jusqu'à la bonne ligne :
			ligne[0].scrollIntoView({block: "center", inline: "start"});
		}
	});

	// Create a new catalog:
	const catalog = aladin.createCatalog({
		name      : "Objects",
		sourceSize: 12,
		color     : '#ff0000',
		shape     : drawFunction
	});
	aladin.addCatalog(catalog);

	// Insert sources inside this new catalog:
	if (withInit)
	{
		refreshAladinLiteCatalog(catalog, 1000);
	}

	return catalog;
}

/**
 * Met a jour le catalogue donne. Le catalogue est entierement efface puis les
 * les 'limit' premiers objets du DataTable y sont inseres.
 * 
 * @param {Object} catalog Le catalogue a mettre a jour.
 * @param {Number} [limit] Nombre maximum d'objets a inserer dans le catalogue.
 *                         Si omis, la limite sera le nombre d'objets visibles
 *                         dans le DataTable.
 */
function refreshAladinLiteCatalog(catalog, limit){
	// empty the catalog:
	catalog.removeAll();

	// get the first objects:
	const coords = getCoord(limit);
			
	// collect objects into an array of AladinSource with more metadata:
	mySrc = [];
	coords.forEach(coord => {
		const pm = getPM(coord.idx);
		mySrc.push(
			aladin.createSource(
				coord.ra,
				coord.dec,
				{
					idx         : coord.idx,
					otype       : getOtype(coord.idx),
					pmra        : (pm ? pm[0] : undefined),
					pmdec       : (pm ? pm[1] : undefined),
					galdim_angle: getGaldimAngle(coord.idx)
				}
			)
		);
	});

	// add all these AladinSources to the Aladin catalog:
	catalog.addSources(mySrc);
}


/**
 * Fonction qui ouvre la fenetre contenant AladinLite si par defaut
 * Lancé au document.ready
 */

function initAladinLiteDialog(s, center, zoom) {
	      if(isAladinLiteWanted()) {
		openAladinLiteDialog(s, center, zoom,true);
		}
}

/**
 * Fonction qui se déplace dans AladinLite et affiche le rayon correspondant au zoom
*/
function moveAladinLitePreview(center, zoom) {
	if (aladin) {
	  aladin.gotoObject(center);
	  aladin.setZoom(zoom*2.5);
	  coo = new Coo();
          if (!coo.parse(center)) return;
	  if (searchCircle!=null) {
		  searchCircle.hide(); 
		  searchCircle.setRadius(zoom);
		  searchCircle.setCenter([coo.lon, coo.lat]);
	  	  searchCircle.show();	
	  }
	  if (simbadCatalogStick) showSimbadCone();
	}
}

/**
 * Fonction qui ouvre la fenetre contenant AladinLite, se déplace et affiche un cercle
 */
function openAladinLitePreview(s, center, zoom) {
	zoom=readRadius();
	openAladinLiteDialog(s, center, zoom*2.5, false);
	simbadCatalogStick=true;
	// certainement que l'objet aladin ne sera pas fini d'etre créé au moment du 1er appel, alors on n'aura pas le rayon, tant pis
	//moveAladinLitePreview(center, zoom);
	changeALDebounced();
	$('#aladin-lite-div').css({'width':'90%','height':'90%','min-width':'250px','min-height':'250px'});
}

/**
 * Fonction qui ouvre la fenetre contenant AladinLite
 */
function openAladinLiteDialog(s, center, zoom, inList) {
	if (aladin==null) {
		if (inList) {
			myAladinLiteList(s, center, zoom);
			$('#aladin-lite-div').css({'height':'90%','min-width':'250px','min-height':'250px'});
		}
		else {
			myAladinLite(s, center, zoom, false);
			// set size of aladin-lite div
			$('#aladin-lite-div').css({'height':'90%','min-width':'250px','min-height':'250px'});
			
		}
		if(useStorage) { 
			var old = localStorage.getItem("ALposition");
			var left = 0;
			if (old) { 
			left = old; 
			$("#aladinLiteDialog").dialog({resizable:true, position:{my:"left top", at:"left+"+left+"px top", of: window}, height:'auto',width:250});
			}
			else { 
				if (inList)
					$("#aladinLiteDialog").dialog({resizable:true, position:{my:"center", at:"center center", of: "#datatable"}, height:'auto',width:250});
				else	$("#aladinLiteDialog").dialog({resizable:true, position:{my:"center", at:"center center", of: window}, height:'auto',width:250});
			}
			if (aladin) aladin.view.requestRedraw();
		}
		else { 
			$("#aladinLiteDialog").dialog({resizable:true, position:{my:"center", at:"center center", of: "#datatable"}, height:'auto',width:250});
		}
		$("#aladinLiteDialog").dialog({
			drag: function(ev, ui) {
				if (useStorage)	localStorage.setItem("ALposition", ui.position.left);
			}
		});

		$("#aladinLiteDialog").parent().css({position:"fixed", top:0, right:0});
	}
	$("#aladinLiteDialog").dialog( "open");
	localStorage.setItem(keyOpen, "true");
}

/**
 * Fonction qui retourne l'index d'une colonne en fonction de son nom
 * @param  {String} name Le nom de la colonne à rechercher
 * @return {Integer}      L'index de la colonne
 */
function getIndexOf(name) {
    var colonnes = $("#datatable th").map(function(){
      return $(this).text();
    }).get();
    for (var i = 0; i < colonnes.length; i++) {
      if(colonnes[i].indexOf(name) > -1) {
        return i;
      }
    }
    return -1;
}


/**
 * Fonction qui place l'icone de prévisualisation dans
 * AladinLite
 */
function setPreviewIcon() {
	var index_col_id = getIndexOf("Identifier");
	var balise_img = "<img class='imgpreview' title='Preview in AladinLite' src='../icons/AladinIconSS.gif' alt=''/>";

	// pour chaque ligne :
	$("#datatable tr").each(function() {

		// ajout de l'image oeil :
		$(this).find('td').eq(index_col_id).append(balise_img);

		// ajout de l'evenement hover :
		$(this).find('td').eq(index_col_id).hover(function() {
				if (localStorage.getItem(keyOpen)== "true") {
					$(this).find("img[class='imgpreview']").css("visibility", "visible");
				}
			},
			function() {
				if(!$(this).parent().hasClass("selectedRow")) {
					$(this).find("img[class='imgpreview']").css("visibility", "hidden");
				}
			}
		);
	});
}


/**
 * Fonction qui surligne dans AladinLite la ligne sur laquelle
 * pointe la souris
 */
function highlightLineHover() {
	var lignes = $("#datatable tr");
	var curSource;

	// highlight in aladin :
	lignes.mouseenter(function() {
		// get the object index:
		const idx = $(this).attr('data-idx');
		if (idx>=0)
		{
			// unselect the current selected source:
			if (curSource)
				curSource.deselect();

			// search for the corresponding source in AladinLite:
			const source = mySrc.find(src => src.data.idx == idx);

			// when found...
			if (source) {
				// ...select it in AladinLite:
				source.select();
				// ...and remember about this new source selection:
				curSource = source;
			}
		}
	});
	lignes.mouseleave(function() {
		if (curSource) {
			curSource.deselect();
		}
	});
}

/**
 * Fonction qui cherche la ligne du tableau representant l'objet objIndex
 * (cad tr[data-idx=objIndex]).
 * 
 * @param {number} objIndex Indice de l'objet a chercher.
 * 
 * @returns {JQueryObject} La ligne du tableau correspondante,
 *                         ou NULL si pas trouvee.
 */
function getTableRow(objIndex){
	if (objIndex === null || typeof(objIndex) === 'undefined' || objIndex < 0)
		return null;
	else
		return document.querySelector(`#datatable tbody tr[data-idx="${objIndex}"]`);
}

/**
 * Fonction qui cherche la cellule du tableau representant l'objet objIndex
 * (cad tr[data-idx=objIndex]) et la colonne colIndex.
 * 
 * @param {number} objIndex Indice de l'objet a chercher.
 * @param {number} colIndex Indice de la colonne a recuperer.
 * 
 * @returns {JQueryObject} La cellule du tableau correspondante,
 *                         ou NULL si pas trouvee.
 */
function getTableCell(objIndex, colIndex){
	if (objIndex === null || typeof(objIndex) === 'undefined' || objIndex < 0
	    || colIndex === null || typeof(colIndex) === 'undefined' || colIndex < 0)
		return null;
	else
		return document.querySelector(`#datatable tbody tr[data-idx="${objIndex}"] td:nth-child(${colIndex+1})`);
}

/**
 * Fonction qui construit et retourne un tableau contenant
 * les coordonnées de chaque enregistrement Simbad
 * @param {number} limit  Nombre maximum de lignes a recuperer.
 * @return {{idx: number, ra: number, dec: number}[]} Ensemble des coordonnees
 */
 function getCoord(limit) {
	// If no limit, use the datatable's page size:
	if(typeof(limit) == "undefined") {
	  limit = (mytable && mytable.page && mytable.page.len() > 0 ? mytable.page.len() : -1);
	}
  
	// Identify RA and DEC columns in the table:
	const index_col_ra = getIndexOf("RA");
	const index_col_dec = getIndexOf("DEC");
  
	let coord = [];
	let i = 0;
	$("#datatable tbody").find("tr").each(function() {
	  // fetch the RA and DEC raw values:
	  const ra_str = $(this).find("td").eq(index_col_ra).html();
	  const dec_str = $(this).find("td").eq(index_col_dec).html();
  
	  // convert the raw coordinates in numeric ones:
	  const coo = new Coo();
	  coo.parse(ra_str + " " + dec_str.substr(1));
  
	  // push the new coordinates:
	  coord[coord.length] = {
		"idx": $(this).attr('data-idx'),
		"ra" : coo.lon,
		"dec": coo.lat
	  };
  
	  // prepare for the next row:
	  i++;
  
	  // stop after 'limit' items:
	  if(i === limit) {
		  return coord;
	  }
	});
  
	return coord;
  }


/**
 * Fonction qui extrait le type de l'objet objIndex dans le tableau.
 * 
 * @param objIndex Numéro d'un objet du tableau de résultats affiché.
 * 
 * @return {string} Type d'objet associé,
 *                  ou '?' si aucun trouvé.
 */
function getOtype(objIndex) {
	// Get the column index of the Otype:
	const index_col_otype = getIndexOf("Otype");
	
	// Get the cell containing the Otype at the given row index:
	const col_otype = getTableCell(objIndex, index_col_otype);

	// Return the otype string, if the cell is found:
	return (col_otype ? col_otype.innerText : '?');
  }

/**
 * Fonction qui extrait l'angle d'orientation galactique pour l'objet objIndex
 * dans le tableau.
 * 
 * @param objIndex Numéro d'un objet du tableau de résultats affiché.
 * 
 * @return {number} Angle d'orientation galactique associé,
 *                  ou NULL si aucun angle trouvé.
 */
 function getGaldimAngle(objIndex) {
	// Get the column index of the Galdim:
	const index_col_galdim = getIndexOf("Angular size");
	
	// Get the cell containing the Galdim at the given row index:
	const col_galdim = getTableCell(objIndex, index_col_galdim);

	// Extract and return only the angle:
	if (col_galdim){
		let parts = col_galdim.innerText.split(' ');
		if (parts && parts.length >= 3){
			const angle = Number.parseFloat(parts[2]);
			if (Number.isFinite(angle))
				return angle;
		}
	}

	// If no Galdim info, nothing to return:
	return null;
  }

/**
 * Fonction qui extrait le mouvement propre de l'objet objIndex dans le tableau.
 * 
 * @param objIndex Numéro d'un objet dans le tableau de résultats affiché.
 * 
 * @return {number} Mouvement propre associé,
 *                  ou NULL si pas trouvé.
 */
function getPM(objIndex) {
	// Get the column index of the Proper Motion:
	const index_col_pm = getIndexOf("Proper motions");
	
	// Get the cell containing the Proper Motion at the given row index:
	const col_pm = getTableCell(objIndex, index_col_pm+1);

	// Extract and return the Proper Motion:
	if (col_pm){
		let parts = col_pm.innerText.split(' ');
		if (parts && parts.length === 2){
			return [ Number.parseFloat(parts[0]), Number.parseFloat(parts[1]) ];
		}
	}

	// If no PM info, nothing to return:
	return null;
}

/**
 * Fonction qui permet de faire le lien entre les enregistrements de Simbad
 * et AladinLite
 */
function linkSimbadWithAladinLite() {
	highlightLineHover();
	setPreviewIcon();

	var images = $(".imgpreview");

	// on update AladinLite lors d'un clic :
	images.click(function() {
		if( !$("#aladinLiteDialog").dialog("isOpen")) {
			openAladinLiteDialog();
			//localStorage.setItem("wantedAL", "true");
		}

		var ligne = $(this).parent().parent();

		selectThisLineInAladin(ligne);
		selectLine(ligne);
	});
}

function  initSurveys(ra, dec) {
	// http://alasky.u-strasbg.fr/MocServer/query?RA=10.8&DEC=32.2&SR=1.5&id=*2MASS*color*
	$.ajax({
		url: '//alasky.u-strasbg.fr/MocServer/query',
		type:'GET',
		cache:true,
		data: {RA: ra, DEC: dec, SR:0.5, id:'*/P/*'},
		complete: function(xhr, textStatus) {
			result  = xhr.responseText;
			if (textStatus!="success") {
				console.log("erreur reponse MocServer :"+xhr);
			}
			else {
			// une fois qu'on a le resultat, grise le menu si besoin
				$('#surveysList option').each ( function (index, element) {
					name=this.text;
					if (result.indexOf(getSurveyName(name)) != -1) {
						$(this).removeAttr('disabled');
					}
					else {
						$(this).attr("disabled",'disabled');
					}
				});
			}
		}
	});
}

function  initSurveysOld(ra, dec) {
	// http://alasky.u-strasbg.fr/MocServer/query?RA=10.8&DEC=32.2&SR=1.5&id=*2MASS*color*
	$('#surveysList option').each ( function (index, element) {
		name=this.text;
		var actuel = $(this);
		$.ajax({
			url: '//alasky.u-strasbg.fr/MocServer/query',
			type:'GET',
			cache:true,
			data: {RA: ra, DEC: dec, SR:1.5, id:'CDS/'+getSurveyName(name)},
			complete: function(xhr, textStatus) {
				if (textStatus!="success") {
					console.log("erreur "+xhr);
				}
				else {
				// une fois qu'on a le resultat, grise le menu si besoin
					result  = xhr.responseText;
					if (result.trim().length==0) {
						actuel.attr("disabled",'disabled');
					}
					else {
						console.log(result);
						actuel.removeAttr('disabled');
					}
				}
			}	
		});
	});
}

function showSimbadCone() {
	var pos=aladin.getRaDec();
	var radius=Math.min(aladin.getFov()[0],1);
	if (simbadCatalog) {
                simbadCatalog.removeAll();
                aladin.view.requestRedraw();
                simbadCatalog=null;
		if (!simbadCatalogStick) return;
        }

	simbadCatalog=A.catalogFromURL(
		'https://simbad.cds.unistra.fr/cone?RA='+pos[0]+'&DEC='+pos[1]+'&SR='+radius+'&MAXREC=30', 
		{sourceSize:12, color:'#30a090'}
	); 
	aladin.addCatalog(simbadCatalog);
	aladin.on('objectClicked', function(object) {   
		if (!object) return;
		var objName=object.data.main_id;
		aladin.showPopup(object.ra,object.dec,'',
		'<a href=\'//simbad.u-strasbg.fr/simbad/sim-id?Ident='+encodeURIComponent(objName)+'\'>'+ objName+'</a>');
	});
	aladin.view.requestRedraw();
}

function addSimbadCone() {

$('#aladin-lite-div').append("\
  <div id='simbad-cone-div' title='Show Simbad objects' onclick='showSimbadCone();'\
 style='top:1.5em;color:darkcyan;cursor:pointer;position: absolute;z-index: 20;background: rgba(250,250,250,0.8);border-radius: 4px;'>\
   <img style='vertical-align:middle' alt='Simbad' src='/icons/simbad.png' height='20'/>\
        </div>\
   ");
$('#aladin-lite-div').prev().css('display','none');
}

function selectSurvey(survey) {
	// Pour 2MASS et DSS, ce sont les radio button
	var name = getShortName(survey);
	if (name=="2MASS" || name=="DSS") {
		$('#'+name).attr('checked', 'checked');
		//document.getElementById(name).checked="checked";
	}
	else {
		document.getElementById("others").checked="checked";
		var liste = document.getElementById("surveysList");

		for (var i = 0 ; i < liste.options.length ; i++) {
			if (liste.options[i].text==name)
				liste.options[i].selected="selected";
		}
	}
	
	// prepare le onchange pour les prochaines fois
	$('#2MASS').on('change', function() {
		surveyTo("2MASS");
	});
	$('#DSS').on('change', function() {
		surveyTo("DSS");
	});
	$('#others').on('change', function() {
		surveyToOthers();
	});
	$('#surveysList').on('change', function() {
		surveyToOthers();
		document.getElementById("others").checked="checked";
	});	
}

function gotoAladinLite(s, id, zoom) {
	if (aladin==null) {
		myAladinLite(s,id,zoom)
	}
	aladin.setImageSurvey(s);
	aladin.setZoom(zoom);
	aladin.gotoObject(decodeURIComponent(id));
	document.getElementById(getShortName(decodeURIComponent(s))).checked="checked";
}

function surveyToOthers() {
	var liste = document.getElementById("surveysList");
	var SelValue = liste.options[liste.selectedIndex].text;
	aladin.setImageSurvey(getSurveyName(SelValue));
}

function surveyTo(s) {
	aladin.setImageSurvey(getSurveyName(s));
}

function surveyTo2MASS() {
	aladin.setImageSurvey("P/2MASS/color");
}

function surveyToDSS() {
	aladin.setImageSurvey("P/DSS2/color");
}

function surveyToXMM() {
	aladin.setImageSurvey("P/XMM/PN/color");
}

function surveyToAllWISE() {
	aladin.setImageSurvey("P/allWISE/color");
}

function surveyToSDSS() {
	aladin.setImageSurvey("P/SDSS9/color");
}

var surveyId2Name = {
  "P/2MASS/color": "2MASS",
  "P/DSS2/color": "DSS",
  "P/XMM/PN/color": "XMM",
  "P/SDSS9/color": "SDSS",
  "P/IRIS/color": "IRIS",
  "P/GALEXGR6/AIS/color": "GALEX",
  "P/allWISE/color" : "AllWISE",
  "P/SPITZER/color": "IRAC" 
};

var name2SurveyId = {};
for (var id in surveyId2Name) {
    name2SurveyId[surveyId2Name[id]] = id;
}

function getShortName(surveyId) {
	return surveyId2Name[surveyId];
}

function getSurveyName(shortName) {
	return name2SurveyId[shortName];
}


/* DIVIDED CONTENT RESIZING FUNCTIONS */

function initSplitTableAladin(eltDatatable){
	// Create and insert the divider wrapper:
	const divWrapper = document.createElement('div');
	divWrapper.setAttribute('class', "divided");
	// ...insert it just before the datatable:
	eltDatatable.parentElement.insertBefore(divWrapper, eltDatatable);

	// Wrap the datatable:
	const tableWrapper = document.createElement('div');
	tableWrapper.setAttribute('class', "divided-content");
	tableWrapper.setAttribute('style', "width: 70vw");
	tableWrapper.appendChild(eltDatatable);
	divWrapper.appendChild(tableWrapper);

	// Append the divider:
	const divider = document.createElement('span');
	divider.setAttribute('class', "divider");
	divider.setAttribute('title', "Change the size of the table and the AladinLite views.");
	divider.innerHTML = '<a href="" class="to-right" title="Hide/Show AladinLite"></a><a href="" class="to-left"  title="Hide/Show table"></a>';
	divWrapper.appendChild(divider);

	// Append AladinLite:
	const aladinWrapper = document.createElement('div');
	aladinWrapper.setAttribute('class', "divided-content AladinLiteBox");
	divWrapper.appendChild(aladinWrapper);

	// Create AladinLite itself:
	const aladinLite = document.createElement('div');
	aladinLite.setAttribute('id', "aladin-lite-div");
	aladinWrapper.appendChild(aladinLite);

	// Init events on these divided content:
	initDividedContent();
}

/**
 * Initialize events of the Divider bar (between the Datatable and AladinLite).
 */
function initDividedContent(){

	/** Function to resize both parts of the division. */
	function resizeDividedContent(evt){
		const viewportWidth   = document.documentElement.clientWidth;
		const viewportPercent = 100*evt.pageX/viewportWidth;
		const width_delta = (evt.pageX - original_mouse_x);
		const aladinWidth = original_width[1] - width_delta;
		/* If AladinLite is becoming too small, hide it completely and show only
		 * the Datatable:
		 * (this aims to fix a bug of AladinLite when its width is 0) */
		if (aladinWidth <= 50){
			// Datatable:
			elements[0].style.width   = '100vw';
			// AladinLite:
			elements[1].style.display = 'none';
			elements[1].style.width   = '50px';
			// Remember AladinLite is hidden for the next times:
			if (useStorage) localStorage.setItem(keyOpen, "false");
		}
		/* Otherwise, just apply the computed relative width: */
		else{
			// Datatable:
			elements[0].style.width = `${viewportPercent}vw`;
			// AladinLite:
			elements[1].style.display = null;
			elements[1].style.width = `${100-viewportPercent}vw`;
			// Remember AladinLite is visible for the next times:
			if (useStorage) localStorage.setItem(keyOpen, "true");
		}
	}

	function stopResizeDividedContent(evt){
	  window.removeEventListener('mousemove', resizeDividedContent);
	  window.removeEventListener('mouseup', stopResizeDividedContent);
	  window.dispatchEvent(new Event('resize'));
	}

	document.querySelector('.divided .divider')
	        .addEventListener('mousedown', evt=>{
	  evt.preventDefault()
	  elements = [...evt.target.parentElement.querySelectorAll('.divided-content')];
	  original_width = elements.map(elt => parseFloat(getComputedStyle(elt, null).getPropertyValue('width').replace('px', '')));
	  original_x = elements.map(elt => elt.getBoundingClientRect().left);
	  original_mouse_x = evt.pageX;
	  window.addEventListener('mousemove', resizeDividedContent);
	  window.addEventListener('mouseup', stopResizeDividedContent);
	});

	document.querySelector('.divided .divider .to-left')
			.addEventListener('click', evt=>{
	  evt.stopPropagation();
	  evt.stopImmediatePropagation();
	  evt.preventDefault();
	  const contents = evt.target.parentElement.parentElement.parentElement.querySelectorAll('.divided-content');
	  contents[0].style.width   = '0px';
	  contents[1].style.display = null;
	  contents[1].style.width   = '100vw';
	  window.dispatchEvent(new Event('resize'));
	  // Remember AladinLite is visible for the next times:
	  if (useStorage) localStorage.setItem(keyOpen, "true");
	});
	document.querySelector('.divided .divider .to-right')
			.addEventListener('click', evt=>{
	  evt.stopPropagation();
	  evt.stopImmediatePropagation();
	  evt.preventDefault();
	  const contents = evt.target.parentElement.parentElement.parentElement.querySelectorAll('.divided-content');
	  contents[0].style.width   = '100vw';
	  contents[1].style.display = 'none';
	  contents[1].style.width   = '50px';
	  // Remember AladinLite is hidden for the next times:
	  if (useStorage) localStorage.setItem(keyOpen, "false");
	});

	// Hide AladinLite if hidden the last time the page was loaded:
	if (!isAladinLiteWanted()){
		document.querySelector('.divided .divider .to-right').click();
	}
}


/* TOOL FUNCTIONS (for broader browser support) */

function hasClass(el, className) {
	return el.classList ? el.classList.contains(className) : new RegExp('\\b'+ className+'\\b').test(el.className);
}

function addClass(el, className) {
	if (el.classList) el.classList.add(className);
	else if (!hasClass(el, className)) el.className += ' ' + className;
}

function removeClass(el, className) {
	if (el.classList) el.classList.remove(className);
	else el.className = el.className.replace(new RegExp('\\b'+ className+'\\b', 'g'), '');
}

function addEvent(el, type, handler) {
    if (el.attachEvent) el.attachEvent('on'+type, handler); else el.addEventListener(type, handler);
}
function removeEvent(el, type, handler) {
    if (el.detachEvent) el.detachEvent('on'+type, handler); else el.removeEventListener(type, handler);
}
