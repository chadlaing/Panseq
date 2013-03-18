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
  $("#query-source-list, #query-target-list").selectable();
  $("#query-add-button").click(queryAdd);
  $("#query-remove-button").click(queryRemove);
  $("#query-remove-all-button").click(queryRemoveAll);

  $("#reference-source-list, #reference-target-list").selectable();
  $("#reference-add-button").click(referenceAdd);
  $("#reference-remove-button").click(referenceRemove);
  $("#reference-remove-all-button").click(referenceRemoveAll);

  //hover events for list items
 $("#query-source-list, #query-target-list").selectable({
  unselected: function(){
    $(":not(.ui-selected)", this).each(function(){
      $(this).removeClass('ui-state-highlight');
    });
  },
  selected: function(){
    $(".ui-selected", this).each(function(){
      $(this).addClass('ui-state-highlight');
    });
  }
});

$(".listContainer li").hover(
  function () {
    $(this).addClass('ui-state-hover');
  }, 
  function () {
    $(this).removeClass('ui-state-hover');
  }
);

}); //end document ready


  // //add and add-all buttons
  function queryAdd() {
    queryTransfer($("#query-source-list li.ui-selected"));
  }

  function queryTransfer(listItems) {
    listItems.fadeOut(0,function() {
      $(this)
        .removeClass("ui-selected")
        .clone()
        .appendTo("#query-target-list")
        .fadeIn(0)
        .data("index", $("#query-source-list li").index($(this)));
    });
  }

  //remove button
  function queryRemove() {
  // $("#target-list li.ui-selected").fadeOut(function() {
    $("#query-target-list li.ui-selected").fadeOut(0,function() {
    $("#query-source-list li")
      .eq($(this).data("index"))
      .removeClass("ui-selected")
      .fadeIn(0);    
     $(this).remove();
    });
  }


  //remove all
  function queryRemoveAll() {
  $("#query-target-list li").fadeOut(0)
    .promise().done(function() {
      $("#query-target-list li").remove();
      $("#query-source-list li:hidden").fadeIn();
    });
  }

  // //add and add-all buttons
  function referenceAdd() {
    referenceTransfer($("#reference-source-list li.ui-selected"));
  }

  function referenceTransfer(listItems) {
    listItems.fadeOut(0,function() {
      $(this)
        .removeClass("ui-selected")
        .clone()
        .appendTo("#reference-target-list")
        .fadeIn(0)
        .data("index", $("#reference-source-list li").index($(this)));
    });
  }

  //remove button
  function referenceRemove() {
  // $("#target-list li.ui-selected").fadeOut(function() {
    $("#reference-target-list li.ui-selected").fadeOut(0,function() {
    $("#reference-source-list li")
      .eq($(this).data("index"))
      .removeClass("ui-selected")
      .fadeIn(0);    
     $(this).remove();
    });
  }

  //remove all
  function referenceRemoveAll() {
  $("#reference-target-list li").fadeOut(0)
    .promise().done(function() {
      $("#reference-target-list li").remove();
      $("#reference-source-list li:hidden").fadeIn();
    });
  }


