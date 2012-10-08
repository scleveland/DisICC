class DisorderValuesController < ApplicationController
  # GET /disorder_values
  # GET /disorder_values.xml
  def index
    @disorder_values = DisorderValue.all

    respond_to do |format|
      format.html # index.html.erb
      format.xml  { render :xml => @disorder_values }
    end
  end

  # GET /disorder_values/1
  # GET /disorder_values/1.xml
  def show
    @disorder_value = DisorderValue.find(params[:id])

    respond_to do |format|
      format.html # show.html.erb
      format.xml  { render :xml => @disorder_value }
    end
  end

  # GET /disorder_values/new
  # GET /disorder_values/new.xml
  def new
    @disorder_value = DisorderValue.new

    respond_to do |format|
      format.html # new.html.erb
      format.xml  { render :xml => @disorder_value }
    end
  end

  # GET /disorder_values/1/edit
  def edit
    @disorder_value = DisorderValue.find(params[:id])
  end

  # POST /disorder_values
  # POST /disorder_values.xml
  def create
    @disorder_value = DisorderValue.new(params[:disorder_value])

    respond_to do |format|
      if @disorder_value.save
        format.html { redirect_to(@disorder_value, :notice => 'Disorder value was successfully created.') }
        format.xml  { render :xml => @disorder_value, :status => :created, :location => @disorder_value }
      else
        format.html { render :action => "new" }
        format.xml  { render :xml => @disorder_value.errors, :status => :unprocessable_entity }
      end
    end
  end

  # PUT /disorder_values/1
  # PUT /disorder_values/1.xml
  def update
    @disorder_value = DisorderValue.find(params[:id])

    respond_to do |format|
      if @disorder_value.update_attributes(params[:disorder_value])
        format.html { redirect_to(@disorder_value, :notice => 'Disorder value was successfully updated.') }
        format.xml  { head :ok }
      else
        format.html { render :action => "edit" }
        format.xml  { render :xml => @disorder_value.errors, :status => :unprocessable_entity }
      end
    end
  end

  # DELETE /disorder_values/1
  # DELETE /disorder_values/1.xml
  def destroy
    @disorder_value = DisorderValue.find(params[:id])
    @disorder_value.destroy

    respond_to do |format|
      format.html { redirect_to(disorder_values_url) }
      format.xml  { head :ok }
    end
  end
end
