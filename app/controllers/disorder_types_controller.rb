class DisorderTypesController < ApplicationController
  # GET /disorder_types
  # GET /disorder_types.xml
  def index
    @disorder_types = DisorderType.all

    respond_to do |format|
      format.html # index.html.erb
      format.xml  { render :xml => @disorder_types }
    end
  end

  # GET /disorder_types/1
  # GET /disorder_types/1.xml
  def show
    @disorder_type = DisorderType.find(params[:id])

    respond_to do |format|
      format.html # show.html.erb
      format.xml  { render :xml => @disorder_type }
    end
  end

  # GET /disorder_types/new
  # GET /disorder_types/new.xml
  def new
    @disorder_type = DisorderType.new

    respond_to do |format|
      format.html # new.html.erb
      format.xml  { render :xml => @disorder_type }
    end
  end

  # GET /disorder_types/1/edit
  def edit
    @disorder_type = DisorderType.find(params[:id])
  end

  # POST /disorder_types
  # POST /disorder_types.xml
  def create
    @disorder_type = DisorderType.new(params[:disorder_type])

    respond_to do |format|
      if @disorder_type.save
        format.html { redirect_to(@disorder_type, :notice => 'Disorder type was successfully created.') }
        format.xml  { render :xml => @disorder_type, :status => :created, :location => @disorder_type }
      else
        format.html { render :action => "new" }
        format.xml  { render :xml => @disorder_type.errors, :status => :unprocessable_entity }
      end
    end
  end

  # PUT /disorder_types/1
  # PUT /disorder_types/1.xml
  def update
    @disorder_type = DisorderType.find(params[:id])

    respond_to do |format|
      if @disorder_type.update_attributes(params[:disorder_type])
        format.html { redirect_to(@disorder_type, :notice => 'Disorder type was successfully updated.') }
        format.xml  { head :ok }
      else
        format.html { render :action => "edit" }
        format.xml  { render :xml => @disorder_type.errors, :status => :unprocessable_entity }
      end
    end
  end

  # DELETE /disorder_types/1
  # DELETE /disorder_types/1.xml
  def destroy
    @disorder_type = DisorderType.find(params[:id])
    @disorder_type.destroy

    respond_to do |format|
      format.html { redirect_to(disorder_types_url) }
      format.xml  { head :ok }
    end
  end
end
