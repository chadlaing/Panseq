$(document).ready(function() {
   
   // the hoverable tabs with jQueryUI
   $(function() {
    $( "#tabs" ).tabs({
      event: "click"
    });
  });


   // submit button with jQueryUI
   $(function() {
    $( "input[type=submit], a, button" )
      .button()
      .click(function( event ) {
        // event.preventDefault();
      });
  });





   
});
