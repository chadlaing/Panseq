#!/usr/bin/env perl
package Interface::Scripts::Panseq;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";
use parent 'Interface::Scripts::Panseq_Super';
use Interface::Scripts::Common;


sub setup{
	my $self=shift;

	$self->start_mode('display');
	$self->run_modes(
		'display'=>'display'
	);
}

sub display{
	my $self=shift;

	my $common = Interface::Scripts::Common->new();

	my $display = $common->_getPageTop();
	$display .= $self->_getPanseqPage();
	$display .= $common->_getPageBottom();
	
	return $display;
}


sub _getPanseqPage{
	my $self =shift;
	my $html = q(

<div class="nav">
	<ul>
		<li><a href="/home/">Home</a></li>
		<li><a href="/panseq/">Panseq</a></li>
		<li><a href="/tutorial/">Tutorial</a></li>
	</ul>
</div>

<div class="columnContainer">
	<a href='#userPan'>
	<div>
		<p>Pan-genome Analyses</p>
	</div>
	</a>


	<a href="#userNovel">
	<div>
		<p>Novel Region Finder</p>
	</div>
	</a>

	<a href="#userLoci">
	<div>
		<p>Loci Selector</p>
	</div>
	</a>
</div>  <!-- End of columnContainer div -->

<div class="userOptions">
	<div id="userPan">
		<form 	method="post",
				enctype="application/x-www-form-urlencoded",
				action="/assays/r"
				class="assayForm"
				>
		 		<fieldset>
		  			<legend>Program Options</legend>
		  			<p>Form 1</p>					 
		 		</fieldset>
		</form>
		<button>Submit</button>
	</div>


	<div id="userNovel">
		<form 	method="post",
				enctype="application/x-www-form-urlencoded",
				action="/assays/r"
				class="assayForm"
				>
		 		<fieldset>
		  			<legend>Program Options</legend>
		  			<p>Form 2</p>

		  			<div id="source-container">
					  <h2 id="source-title">Widgets</h2>
					  <ul id="source-list">
					    <li>Accordion</li>
					    <li>Autocomplete</li>
					    <li>Button</li>
					    <li>Datepicker</li>
					    <li>Dialog</li>
					    <li>Progressbar</li>
					    <li>Slider</li>
					    <li>Tabs</li>
					  </ul>
					</div>

					<div id="target-container">
				  <h2 id="target-title">Widgets You Used</h2>
				  <ul id="target-list">
				  </ul>

				  <div id="transfer-buttons-container">
				  <div id="transfer-buttons">
				    <button id="add-button">Add &rarr;</button>
				    <button id="add-all-button">Add All *&rarr;</button>
				    <button id="remove-button">&larr; Remove</button>
				    <button id="remove-all-button">&larr;* Remove All</button>           
				  </div>
				</div>
	
	<!-- 			$(document).ready(function() {
			    $("#source-list, #target-list").selectable();
			 
			    $("#add-button").click(add);
			    $("#add-all-button").click(addAll);
			    $("#remove-button").click(remove);
			    $("#remove-all-button").click(removeAll);
			 
			    addHiglightPlugin();
			  }); -->





</div>

		 		</fieldset>
		</form>
		<button>Submit</button>
	</div>


	<div id="userLoci">
		<form 	method="post",
				enctype="application/x-www-form-urlencoded",
				action="/assays/r"
				class="assayForm"
				>
		 		<fieldset>
		  			<legend>Program Options</legend>
		  			<p>Form 3</p>					 
		 		</fieldset>
		</form>
		<button>Submit</button>
	</div>

</div>  <!-- End of userOptions div -->



);
	return $html;
}



1;


