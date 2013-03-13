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
        event.preventDefault();
      });
  });

   // init for the list transfer
  $("#source-list, #target-list").selectable();

  $("#add-button").click(add);
  $("#add-all-button").click(addAll);
  $("#remove-button").click(remove);
  $("#remove-all-button").click(removeAll);
  addHiglightPlugin();

}); //end document ready


  // //add and add-all buttons
  function add() {
    transfer($("#source-list li.ui-selected"));
  }

  function addAll() {
    transfer($("#source-list li:visible"));
  }

  function transfer(listItems) {
    listItems.fadeOut(0,function() {
      $(this)
        .removeClass("ui-selected")
        .clone()
        .appendTo("#target-list")
        .fadeIn(0)
        .data("index", $("#source-list li").index($(this)))
        .highlight();
    });
  }

  //remove button
  function remove() {
  // $("#target-list li.ui-selected").fadeOut(function() {
    $("#target-list li.ui-selected").fadeOut(0,function() {
    $("#source-list li")
      .eq($(this).data("index"))
      .removeClass("ui-selected")
      .fadeIn(0)    
      .highlight();
     $(this).remove();
    });
  }


  //remove all
  function removeAll() {
  $("#target-list li").fadeOut(0)
    .promise().done(function() {
      $("#target-list li").remove();
      $("#source-list li:hidden").fadeIn().highlight();
    });
  }

  //highlight plugin
  function addHiglightPlugin() {
  $.fn.highlight = function() {
    return this
      .addClass("li-transfer-highlight")
      .removeClass("li-transfer-highlight", 400);
    }
  }

