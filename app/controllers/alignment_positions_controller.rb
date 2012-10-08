class AlignmentPositionsController < ApplicationController
  # GET /alignment_positions
  # GET /alignment_positions.xml
  def index
    @alignment_positions = AlignmentPosition.all

    respond_to do |format|
      format.html # index.html.erb
      format.xml  { render :xml => @alignment_positions }
    end
  end

  # GET /alignment_positions/1
  # GET /alignment_positions/1.xml
  def show
    @alignment_position = AlignmentPosition.find(params[:id])

    respond_to do |format|
      format.html # show.html.erb
      format.xml  { render :xml => @alignment_position }
    end
  end

  # GET /alignment_positions/new
  # GET /alignment_positions/new.xml
  def new
    @alignment_position = AlignmentPosition.new

    respond_to do |format|
      format.html # new.html.erb
      format.xml  { render :xml => @alignment_position }
    end
  end

  # GET /alignment_positions/1/edit
  def edit
    @alignment_position = AlignmentPosition.find(params[:id])
  end

  # POST /alignment_positions
  # POST /alignment_positions.xml
  def create
    @alignment_position = AlignmentPosition.new(params[:alignment_position])

    respond_to do |format|
      if @alignment_position.save
        format.html { redirect_to(@alignment_position, :notice => 'Alignment position was successfully created.') }
        format.xml  { render :xml => @alignment_position, :status => :created, :location => @alignment_position }
      else
        format.html { render :action => "new" }
        format.xml  { render :xml => @alignment_position.errors, :status => :unprocessable_entity }
      end
    end
  end

  # PUT /alignment_positions/1
  # PUT /alignment_positions/1.xml
  def update
    @alignment_position = AlignmentPosition.find(params[:id])

    respond_to do |format|
      if @alignment_position.update_attributes(params[:alignment_position])
        format.html { redirect_to(@alignment_position, :notice => 'Alignment position was successfully updated.') }
        format.xml  { head :ok }
      else
        format.html { render :action => "edit" }
        format.xml  { render :xml => @alignment_position.errors, :status => :unprocessable_entity }
      end
    end
  end

  # DELETE /alignment_positions/1
  # DELETE /alignment_positions/1.xml
  def destroy
    @alignment_position = AlignmentPosition.find(params[:id])
    @alignment_position.destroy

    respond_to do |format|
      format.html { redirect_to(alignment_positions_url) }
      format.xml  { head :ok }
    end
  end
end
