//// SAMP integration
// dependance avec samp.js
//cachedScript('http://astrojs.github.io/sampjs/examples/lib/samp.js').done( function() {
$( document ).ready(function() {
cachedScript('/Simlib/samp.js').done( function() {

     // creation de l'object de connexion
        var connector = new samp.Connector("Simbad");
    var hubAvailable = false;
    // Toutes les 2 sec, cherche un hub
    //onload = function() {
    //  connector.onHubAvailability(function(hubRunning) {hubAvailable = hubRunning;}, 2000);
    //};
    // A la deconnection de la page, deconnecte le hub
    onunload = function() {
      connector.unregister();
    };
    // Action du click sur le bouton
    $('[name=sendBySAMP]').click(function() {
	// cherche un hub
        //connector.onHubAvailability(function(hubRunning) {hubAvailable = hubRunning;}, 200);
	samp.ping(function(hubRunning) {hubAvailable = hubRunning;});
	setTimeout(function () {
        if ( ! hubAvailable) {
            alert('No hub available. Launch Topcat or Aladin ...');
            return;
        }
	// Se connecte au hub et envoie le message
        connector.runWithConnection(function(connection) {
	    // 1er message pour avoir les points simbad (URL de la page en sortie votable)
            var msg = '';
	    if (document.URL.indexOf("?")==-1) { // cas ou pas d'interro Votable possible, par ex erreur
	        msg = new samp.Message("script.aladin.send", {"script": "get Simbad"});
	    }
	    else {
	        msg = new samp.Message("table.load.votable", {"url": window.document.URL + "&output.format=VOTable", "name": document.title});
	    }
	    var center = $('#Coord').val();
	    if (typeof center === 'undefined') {
		center = $('#ident').val();
		if (typeof center === 'undefined') {
			center = $("td.lon")[0].innerHTML+$("td.lat")[0].innerHTML;
//		    center = "";
	    	}
	    }
            var radius = $('#radius').val();
            if (typeof radius !== 'undefined') {
                radius = radius*2.5 + $('#radius\\.unit').val();
                center = center + "; zoom "+ radius;
            }
            else {
                center = center + "; zoom 15arcmin";
            }
            // Zoom + centrer : setconf frame=ICRS;00 48 35.40-72 52 56.5;zoom 1.2arcmin
            var msg3 = new samp.Message("script.aladin.send", {"script": "setconf frame=ICRS;"+center});
	    // envoie des messages
            if (document.URL.indexOf("?")==-1) {
	        connection.notifyAll([msg3]);
	  	connection.notifyAll([msg]);
	     }
	    else {
	  	connection.notifyAll([msg]);
	        connection.notifyAll([msg3]);
	     }

//            connection.notifyAll([msg2]);
        });
	},800);
    });

});// end cachedScript samp.js
});// end document.ready


